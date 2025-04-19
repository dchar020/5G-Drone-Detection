

idx = 6;



figure('Renderer', 'painters', 'Position', [10 10 900 600])
subplot(2,3,1);
plot(1:length(y_theo(:,idx)), y_theo(:,idx));

subplot(2,3,4);
plot(1:length(y_exp(:,idx)), y_exp(:,idx));

subplot(2,3,2);
plot(1:length(y_theo(:,idx)), y_theo(:,idx-1));

subplot(2,3,5);
plot(1:length(y_exp(:,idx)), y_exp(:,idx-1));

subplot(2,3,3);
plot(1:length(y_theo(:,idx)), y_theo(:,idx+1));

subplot(2,3,6);
plot(1:length(y_exp(:,idx)), y_exp(:,idx+1));
%}

mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(1:length(y_theo(ridx,:)), y_theo(ridx,:));


%plotResponse(rdresp,y_theo,mfcoeff);

%{
figure;
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(y_theo(ridx,:),PRF,'spectrogram')
%}

M = 80;
L = 16;
lk = 0.7;
x = abs(fftshift(fft(ymf(ridx,1:M), 80)));

figure('Renderer', 'painters', 'Position', [10 10 900 600])
plot(1:length(x), x);



figure('Renderer', 'painters', 'Position', [10 10 900 600])
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y_theo);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range
pspectrum(ymf(ridx,:),PRF,'spectrogram');


%{
, ...
    TimeResolution=M/PRF,OverlapPercent=L/M*100, ...
    Leakage=lk)
%}