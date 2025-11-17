clc; clear; close all;

sim_time = 1e-4;
T = 1e-7;
N_sc = 6;
fs = 2e8;

t = 0:1/fs:sim_time;
t_symbol = 0:1/fs:T;
N = length(t);
cw_waveform = zeros(N_sc, N);
N_fft = N;

for k = 0:N_sc-1
    cw_waveform(k+1,1:length(t_symbol)) = exp(1j*2*pi*k*t_symbol/T);
    cw_waveform(k+1,length(t_symbol):length(t)) = 0;
end

CW_FFT = zeros(N_sc, N_fft);

for k=1:N_sc
    CW_FFT(k,:) = abs(fftshift(fft(cw_waveform(k,:).', N_fft)));
    CW_FFT(k,:) = CW_FFT(k,:)/max(CW_FFT(k,:));
end

fs_range = (-fs/2:fs/N_fft:fs/2-1)/1e6;
figure('Renderer', 'painters', 'Position', [10 10 900 600])
hold on
plot(fs_range, CW_FFT(1,:), 'r');
plot(fs_range, CW_FFT(2,:), 'g');
plot(fs_range, CW_FFT(3,:), 'b');
plot(fs_range, CW_FFT(4,:), 'r');
plot(fs_range, CW_FFT(5,:), 'g');
plot(fs_range, CW_FFT(6,:), 'b');
xlabel('Frequency (MHz)'); ylabel('Normalized Magnitude');
title('Baseband Frequency Spectrum');
hold off



