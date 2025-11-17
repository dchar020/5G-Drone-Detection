

close all;


figure;
plot((1:length(cw_waveform)), cw_waveform);
title('CW 5G OFDM Transmit Signal Pulse')
xlabel('Discrete-time Sample Number'); ylabel('Amplitude')

num_steps = N - length(t_symbol);
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
y = sum(abs(ymf),2);
r = c/(2*BW);
ax = (0:num_steps)*r;
y_ss = (y((length(y)-num_steps):length(y)));
th = num_steps*r+r/2;
a = r/2:r:th;

figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(3,1,1)
plot(ax, y_ss);
grid on
xticks(a);
title('MM Model Range Bin Correlation with Transmit Signal')
xlabel('Range (m)'); ylabel('Correlation')

subplot(3,1,2)
ymf = mf(y_theo_ps);
y = sum(abs(ymf),2);
y_ss = y((length(y)-num_steps):length(y));
plot(ax, y_ss);
grid on
xticks(a);
title('PS Model Range Bin Correlation with Transmit Signal')
xlabel('Range (m)'); ylabel('Correlation')

subplot(3,1,3)
ymf = mf(y_exp);
y = sum(abs(ymf),2);
y_ss = y((length(y)-num_steps):length(y));
plot(ax, y_ss);
grid on
xticks(a);
title('Exp Model Range Bin Correlation with Transmit Signal')
xlabel('Range (m)'); ylabel('Correlation')


%}

idx = 6;

%figure('Renderer', 'painters', 'Position', [10 10 900 600])
figure('Renderer', 'painters', 'Position', [10 10 1000 900])

subplot(3,3,1);
plot(1:length(y_theo(:,idx)), y_theo(:,idx-1));
title('MM model \textbf{r}$_x[:,1]$', Interpreter='latex')
subplot(3,3,2);
plot(1:length(y_theo(:,idx)), y_theo(:,idx));
title('MM model \textbf{r}$_x[:,2]$', Interpreter='latex')
subplot(3,3,3);
plot(1:length(y_theo(:,idx)), y_theo(:,idx+1));
title('MM model \textbf{r}$_x[:,3]$', Interpreter='latex')

subplot(3,3,4);
plot(1:length(y_theo(:,idx)), y_theo_ps(:,idx-1));
title('PS model \textbf{r}$_x[:,1]$', Interpreter='latex')
subplot(3,3,5);
plot(1:length(y_theo(:,idx)), y_theo_ps(:,idx));
title('PS model \textbf{r}$_x[:,2]$', Interpreter='latex')
subplot(3,3,6);
plot(1:length(y_theo(:,idx)), y_theo_ps(:,idx+1));
title('PS model \textbf{r}$_x[:,3]$', Interpreter='latex')

subplot(3,3,7);
plot(1:length(y_exp(:,idx)), y_exp(:,idx-1));
title('Exp model \textbf{r}$_x[:,1]$', Interpreter='latex')
xlabel('Discrete-time Sample Number');
subplot(3,3,8);
plot(1:length(y_exp(:,idx)), y_exp(:,idx));
title('Exp model \textbf{r}$_x[:,2]$', Interpreter='latex')
xlabel('Discrete-time Sample Number');
subplot(3,3,9);
plot(1:length(y_exp(:,idx)), y_exp(:,idx+1));
title('Exp model \textbf{r}$_x[:,3]$', Interpreter='latex')
xlabel('Discrete-time Sample Number');

%}
%{
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
y = sum(abs(ymf),2);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
%}

%{
figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(1:length(y_theo(ridx,:)), y_theo(ridx,:)/abs(max(y_theo(ridx,:))));
xlabel('Discrete Sample Number'); ylabel('Normalized Magnitude');
title('Theoretical MM Time-Domain Radar Signal');
%}
%{
figure('Renderer', 'painters', 'Position', [10 10 900 600])
plotResponse(rdresp,y_theo,mfcoeff);
%}
%{
figure;
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(y_theo(ridx-1,:),PRF,'spectrogram')
%}

%{
M = 80;
L = 16;
lk = 0.7;
x = abs(fftshift(fft(ymf(ridx-1,1:M), 80)));

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(1:length(x), x);
%}


figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(1,3,1)
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(ymf(ridx,:),PRF,'spectrogram');

subplot(1,3,2)
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo_ps);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(ymf(ridx,:),PRF,'spectrogram');

subplot(1,3,3)
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_exp);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(ymf(ridx,:),PRF,'spectrogram');

rdresp  = phased.RangeDopplerResponse('PropagationSpeed',c,'SampleRate',fs,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',128,'DopplerOutput','Speed',...
    'OperatingFrequency',fc);
mf = cw_waveform(length(t_symbol):-1:1);
mfcoeff = conj(mf);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(1,3,1)
plotResponse(rdresp,y_theo,mfcoeff);
%colorbar('Ticks',[-150, -30]);

subplot(1,3,2)
plotResponse(rdresp,y_theo_ps,mfcoeff);

subplot(1,3,3)
plotResponse(rdresp,y_exp,mfcoeff);
%{
, ...
    TimeResolution=M/PRF,OverlapPercent=L/M*100, ...
    Leakage=lk)
%}

%}
figure;
plot(1:16, mfcoeff);

figure; 
temp = abs(y_mf_mm(:,2));
plot(1:30, temp(1:30));

mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo(:,2));
figure;
plot(1:length(ymf), abs(ymf));



figure;
ambgfun(y_mf_mm(:,26),fs,PRF);

