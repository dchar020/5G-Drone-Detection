# Coherent UAV target modelling for OFDM radar

This work allows the predicting of pulsed Orthogonal Frequency Division Multiplexing (OFDM) radar backscatter signal reflected off of dark Uncrewed Aerial Vehicles (UAV)s; ``dark'' meaning the UAV does not actively emit electromagnetic radiation. Modeling radar backscatter from UAV movement in any orientation and position is challenging as the UAV's local positions and micro-velocities must continually be updated as the simulation progresses in time. The proposed mathematical model simplifies the process of constructing these simulations by introducing reference frames anchored to specific points on the UAV and making logical simplifications. The proposed mathematical model allows a simulation designer to construct the radar target in any arbitrary topology to generate micro-Doppler signatures and range-responses with a closed-form solution. The result of the model is compared to standard simulation techniques which have been shown to match real-world data, showing that the model is accurate.

## How to use

The main simulation tool is found in the uav_backscatter_simulator.m file. One may simulate a large variety of scenarios for UAV(s) in range of a monostatic OFDM transeiver. The code is written in a general manner but with the following scenario in-mind: the monostatic transceiver is at [0,0,0] in the coordinate system and one UAV's center-of-mass is located at [100,0,0] meters. Throughout the simulation, the drone only hovers at this fixed location and does not move. The effect of the spinning propellers is shown through the simulation output.

If one would like to calculate the received data matrix from a different scenario, one would have to update the simulation parameters to what is desired (e.g. carrier frequency, number of subcarriers, topology of UAV, number of UAVs, etc). The functions which control how the UAVs move through space must also be updated. These two functions follow:

```matlab

function [Omega,x_p_u,R_p_u,f_rot] = update_prop_params(t_p, i)
    global rpm  Omega_global x_p_u_global R_p_u_global PRF f_rot_global;
    % same rotation frequency for all blades for all sim time.
    f_rot = f_rot_global(i);
    
    % blade angle changes based on the speed of their rotation
    % they all share the same rotation speed so f_rot can be used
    Omega = Omega_global(i) + f_rot/PRF;

    % the position and orientation of the ith propeller wrt the center of 
    % the drone is always the same for this simulation
    x_p_u = x_p_u_global(1:3,i);
    R_p_u = R_p_u_global(1:3,1:3,i);
end

function [x_u_o,x_t_o,v_u_o,v_t_o,R_u_o,f_eps,sigma] = update_params(t_p)
    global R_m RCS_static;
    x_u_o = [R_m;0;0]; 
    x_t_o = [0;0;0];
    v_t_o = [0;0;0];
    v_u_o = [0;0;0];
    R_u_o = eul2rotm([0 0 0]);
    f_eps = 0;
    sigma = RCS_static;
end
```

Each of the function outputs are well defined in both the conference paper (ref) and the thesis (ref). If a user would like to make a change to the UAV's trajectory, for example, they would have to change the behavior of the update_params(t_p) function to return the position of the center-of-mass of the UAV, x_u_o, as well as its velocity, v_u_o, as a function of the time in the simulation (t_p). 

The code offers a way to compare this result with the state-of-the-art simulation point-scatterer methodology through both comparing it with a third-party simulation and with a version which has been implemented by the author. If a similar comparison would be desired by the user, they must implement the code to calculate the positions and velocities of each of the point-scatterers in each of the UAV targets themselves and alter what has been written for the specific scenario. This part may prove challenging which emphasizes the benefit of the novel simulation methodology.

Many figures can be regenerated, such as the micro-Doppler signatures, the range-speed responses, etc. through the other code files: show_fft_slice.m (ref), show_td_sig.m (ref), compare_td_sigs.m (ref). 

# Real-World Data Micro-Doppler Processing

Data has been collected experimentally from spinning a rotor in a bistatic configuration, where the TX and RX channels use separate antennas. The setup is shown in the appendix of the thesis defence presentation slides (ref). The collected data is processed in the processing.m (ref) file. 

## Code explanation

First, one of the radar received data files are loaded into memory. The following code shows which file is which and opens the received data resulting from the propeller spinning at its highest speed.

```matlab
% STATIC
file0 = 'G:\Daniel Test 20251028\File.iq\File_2025-10-29165731.complex.1ch.float32';
% MEDIUM
file1 = 'G:\Daniel Test 20251028\File_001.iq\File_2025-10-29165844.complex.1ch.float32';
% SLOW
file2 = 'G:\Daniel Test 20251028\File_002.iq\File_2025-10-29170243.complex.1ch.float32';
% FAST
file3 = 'G:\Daniel Test 20251028\File_003.iq\File_2025-10-29170840.complex.1ch.float32';

fid = fopen(file3, 'rb');
```

Based on how the data is stored, the real part of the signal is contained in the even index data bins while the imaginary data is in the odd index. These are loaded into the iq data matrix.

```matlab
% I = even index, Q = odd index
N_samples = length(data)/2;
iq = complex(data(1:2:2*N_samples), data(2:2:2*N_samples));
```

Next, based on the speed of the propeller and the sampling frequency during data collection in the experiment, downsampling is performed to increase the frequency resolution of the resulting micro-Doppler signature plot. The _re appendix to the variables indicates the variables have been resampled.

```matlab
downsample_rate = 10;
fs_re = fs/downsample_rate;
iq_re = decimate(iq,downsample_rate);
```

The signal is then altered to remove its mean from each of the data cells. This ensures the signal is better centered around an amplitude of 0 and thus produces more insightful Fourier transforms.

```matlab
iq_re = iq_re - mean(iq_re);
```

Next, we intend to perform time-frequency analysis on the signal through Short Time Fourier Transform (STFT). A large overlap window is chosen to provide smoother transitions between adjacent times. We also allocate twice the window size for the resulting FFT which allows for additional smoothing of the FFT signal. Note that doing this does not increase the precision of the resulting FFT as no additional information is gained. 

```matlab
exponent = 11;
window = hamming(2^exponent);
N_overlap = round(0.95*length(window));
nfft = 2^(exponent+1);

% perform the STFT and plot it

[s,f,ts] = spectrogram(iq_re,window, N_overlap, nfft, fs_re, 'centered');
Sdb = 20*log10(abs(s));
```

The plots are then generated and their frequency and time resolutions are displayed.