clc; clear; close all;

sim_time = 0.1;
fc = 10e3;
BW = 1e3;
fc1 = fc-BW/4;
fc2 = fc+BW/4;
fs = 4*fc;
t = 0:1/fs:sim_time;

x = exp(1j*2*pi*fc1*t) + exp(1j*2*pi*fc2*t);
N = length(t);
fft_spacing = -fs/2+1:fs/N:fs/2;

tx_dft = abs(fftshift(fft(x, N)));
n_tx_dft = tx_dft/max(tx_dft);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(fft_spacing/1e3, n_tx_dft, 'b');
xlabel('Frequency (kHz)'); ylabel('Normalized Magnitude');
xline((fc+BW/2)/1e3); xline((fc-BW/2)/1e3);
title('Frequency Spectrum of Modulated Signal');

%demodulation
x = x.*exp(-1j*2*pi*fc*t);

tx_dft = abs(fftshift(fft(x, N)));
n_tx_dft = tx_dft/max(tx_dft);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(fft_spacing/1e3, n_tx_dft, 'b');
xlabel('Frequency (kHz)'); ylabel('Normalized Magnitude');
xline(+BW/2/1e3); xline(-BW/2/1e3);
title('Frequency Spectrum of Demodulated Signal');
