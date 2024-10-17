% Channel estimation according to the IR1 method


%{
Vraag 1

Modify the m-file IR1.m that conducts a simple experiment to estimate
the impulse response (IR) of the acoustic channel, by literally applying
the definition of the IR (i.e., ‘the response of the system when applying an
impulse at the input’). The m-file should plot a figure with two subplots
containing the estimated IR response (time-domain), and the (magnitude
of the) frequency response. The time-domain scale must be in samples
(‘filter taps’), not seconds, and the frequency scale must be in Hz. The
frequency response is plotted on a dB-scale (on the magnitude axis, not
on the frequency axis). What do you observe?

%}

%% Cleanup
clear; clc; close all

%% Initialize script parameters.
fs = 16000; % Sampling frequency [Hz]
dftsize = 2; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) [samples]
Noverlap = 0; 

%% Create the signal to be played

T = 1/fs;
duration = 2;
n = duration/T; %totaal aantal samples
ni = 0:n-1;
x = 0:T:duration-T; 
y = dirac(x-ni*T);
idx = y == Inf; % find Inf
y(idx) = 1;     % set Inf to finite value
stem(x,y)

sig = [];                       

%% Play and record.
% Call to initparams()

%{
[] = initparams();
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Calculate the IR by trimming the output
h = out();

%% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot();
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
% Magnitude response
subplot(2,1,1)
plot(pow2db());
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')

%% Filtered white noise Vs. recorded white noise.
sig = wgn(); % Generate white noise signal    

% Convolve the IR with the white noise
gen_noise = fftfilt();                  

% Send the white noise across the channel
[] = initparams();
sim('recplay');
rec_noise=simout.signals.values(:,1);

%% Spectrogram and PSD plot
% Spectrogram
spectrogram(); % Generated white noise
spectrogram(); % Recorded white noise

% PSD
[] = spectrogram(); % Generated white noise
[] = spectrogram(); % Recorded white noise

PSD_gen_noise = ; % Generated white noise PSD
PSD_rec_noise = ; % Recorded white noise PSD (ideally only when the signal is active)

figure;
plot(,pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Generated white noise.')
% Output signal
subplot(2,2)
plot(,pow2db());
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Recorded white noise')
%}