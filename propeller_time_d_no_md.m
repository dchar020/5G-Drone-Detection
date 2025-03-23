clear all;

% Tutorial:
% https://www.mathworks.com/help/phased/gs/end-to-end-radar-system.html

% generate the waveform
waveform = phased.RectangularWaveform('PulseWidth',5e-5, ...
    'PRF',5e3,'OutputFormat','Pulses','NumPulses',1);

% display the waveform
%Y = waveform();
%plot(Y)

% define the antenna to use
antenna = phased.IsotropicAntennaElement('FrequencyRange',[1e9 10e9]);

% define the target, propagation characteristics and operating freq
%target = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',[1 1 1 1 1 1 1 1 1 1 1 1], ...
%    'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',4e9);
target = phased.RadarTarget('Model','Nonfluctuating','MeanRCS',[1], ...
    'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',4e9);
bladelen = 4;
scatterpos = [(-2:-2:-6) (2:2:6) 0 0 0 0 0 0;...
              0 0 0 0 0 0 (-2:-2:-6) (2:2:6);...
              zeros(1, 12)];

scatterpos = [702; 0; 0];

% define the position and velocities of the antenna and target
antennaplatform = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
targetplatform = phased.Platform('InitialPosition',[700; 0; 0], ...
    'Velocity',[0;0;0]);

% get the target range and target angles (azim and elev)
[tgtrng,tgtang] = rangeangle(targetplatform.InitialPosition, ...
    antennaplatform.InitialPosition);

% modelling the transmission
numpulses = 1;
SNR = 20;

% SNR here is the required SNR for those params (probability of detection
% and probability of false alarm

% Pt is roughly 45 kW. use 50 kW to be safe. Model the tranmission
% using the phase array object
transmitter = phased.Transmitter('PeakPower',50e3,'Gain',20,'LossFactor',0, ...
    'InUseOutputPort',true,'CoherentOnTransmit',true);

% model the radiator and collector;
% radiator: "A radiator converts signals into radiated wavefields"
% collector is the opposite
radiator = phased.Radiator('Sensor',antenna,...
    'PropagationSpeed',physconst('LightSpeed'),'OperatingFrequency',4e9);
collector = phased.Collector('Sensor',antenna,...
    'PropagationSpeed',physconst('LightSpeed'),'Wavefront','Plane', ...
    'OperatingFrequency',4e9);

% model the receiver:
receiver = phased.ReceiverPreamp('Gain',20,'NoiseFigure',2, ...
    'ReferenceTemperature',290,'SampleRate',1e6, ...
    'EnableInputPort',true,'SeedSource','Property','Seed',1e3);

% model the propagation channel:
channel = phased.FreeSpace(...
    'PropagationSpeed',physconst('LightSpeed'), ...
    'OperatingFrequency',4e9,'TwoWayPropagation',true, ...
    'SampleRate',1e6);

% all the parameters and handles are setup. the following code actions
% them

T = 1/waveform.PRF;
% Get antenna position
txpos = antennaplatform.InitialPosition;
% Allocate array for received echoes
rxsig = zeros(waveform.SampleRate*T,numpulses);

for n = 1:numpulses
    % Update the target position
    [tgtpos,tgtvel] = targetplatform(T);
    % Get the range and angle to the target
    [tgtrng,tgtang] = rangeangle(tgtpos,scatterpos);
    % Generate the pulse
    sig = waveform();
    % Transmit the pulse. Output transmitter status
    [sig,txstatus] = transmitter(sig);
    % Radiate the pulse toward the target
    sig = radiator(sig,tgtang);
    % Propagate the pulse to the target in free space
    %sig = channel(sig,txpos,scatterpos,[0;0;0],zeros(3, 12));
    sig = channel(sig,txpos,scatterpos,[0;0;0],[0;0;0]);
    % Reflect the pulse off the target
    sig = target(sig);
    % Propagate the echo to the antenna in free space
    %sig = channel(sig,scatterpos,txpos,zeros(3, 12),[0;0;0]);
    % Collect the echo from the incident angle at the antenna
    sig = collector(sig,tgtang);
    % Receive the echo at the antenna when not transmitting
    rxsig(:,n) = receiver(sig,~txstatus);
end


%{
res = zeros(200, 1);
for i = 1:numpulses
    %subplot(numpulses, 1, i)
    %plot(1:200, abs(rxsig(:, i)))
    res = res + abs(rxsig(:, i));
end
plot (res)
%}
% pulsint is just adding all the resulting signals for each pulse.
% this instruction reformats rxsig from [200, 10] to [200, 1]
rxsig = pulsint(rxsig,'noncoherent');

% unigrid is uniform grid; a generator for basic expressions
t = unigrid(0,1/receiver.SampleRate,T,'[)');

% range gates are the range bins; there's the same amount of bins as there
% are time intervals for this sampling rate
rangegates = (physconst('LightSpeed')*t)/2;

% plot.
plot(rangegates/1e3,rxsig)
hold on
xlabel('range (km)')
ylabel('Power')
xline(1.4,'r')
hold off


%{
% unigrid is uniform grid; a generator for basic expressions
t = unigrid(0,1/receiver.SampleRate,T,'[)');

% range gates are the range bins; there's the same amount of bins as there
% are time intervals for this sampling rate
rangegates = (physconst('LightSpeed')*t)/2;

% plot.
plot(1:200,rxsig)
hold on
%xline([tgtrng/1e3,tgtrng/1e3],'r')
hold off
%}