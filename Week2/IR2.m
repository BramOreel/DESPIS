% Channel estimation according to the IR2 method

%% Cleanup
clear; clc; close all

%% Initialize script parameters.
fs = 16000; % Sampling frequency [Hz]
dftsize = 100; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) [samples]
Noverlap = dftsize/2; 
channelLength = 200; % Length of impulse response [samples]
delay = 50; % Positive delay safety margin when aligning input and output [samples]

%% Create the signal to be played.
duration = 1; % Duration of the signal in [s]
N =duration*fs
sig = wgn(N,1,0);  

%% Filter signal -> only to be used for exercise 2.3 (To this end, also copy the content
%% of the current file to a new IR_bandstop.m file)
% filt = fir1();
% sig = fftfilt();

%% Play and record.
% Call to initparams()
[simin,nbsecs,fs] = initparams(sig,fs,3);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Calculate the impulse response.
uMatrix = toeplitz(sig'); % Toeplitz matrix


%Find the max of the out signal
c = max(out);
c_60 = 0.6*c;
%Find the first x value that reaches this threshold

for i = 1:size(out,1)
    if(out(i) > c_60)
        x = i;
        break
    end    
end

yOnset = x; % Determine start of recorded signal [samples]

y = out(yOnset-200:yOnset + N-201); % Extract the relevant output signal

h =  lsqr(uMatrix,y)       ; % Estimate impulse response


%% Calculate the IR by trimming the output

%{
c2 = max(out);
x2 = 0;
tol = 10^-9

for i = 1:size(out,1)
    if(abs(out(i) -c2) < tol)
        x2 = i;
        break
    end
    
end

h = out(x2-25:x2+150);
%}
save('channel.mat','h'); % Save impulser response

%% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot(h);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
% Magnitude response
subplot(2,1,2)
plot(pow2db(abs(h).^2));
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')

%{

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