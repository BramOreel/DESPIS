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
delay = 200; % Positive delay safety margin when aligning input and output [samples]

%% Create the signal to be played.
duration = 1; % Duration of the signal in [s]
N = duration*fs;
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

%figure; 
%subplot(3,1,1)
%plot(sig)

%subplot(3,1,2)

%plot(out)

%% Calculate the impulse response.
uMatrix = toeplitz(sig'); % Toeplitz matrix van de ruis


%Find the max of the out signal
c = max(out);
c_60 = 0.6*c;
x = 0;
%Find the first x value that reaches this threshold

for i = 1:size(out,1)
    if(out(i) > c_60)
        x = i;
        break
    end    
end

Onset = x; % Determine start of recorded signal [samples]

%y = out(yOnset-200:yOnset + N-201); % Extract the relevant output signal
y = out(Onset-delay:Onset + N-delay-1);
%subplot(3,1,3)
%plot(y)
%Causaal maken: shiften van out om h causaal te maken:


h =  lsqr(uMatrix,y,0.1)       ; % Estimate impulse response
save('channel.mat','h'); % Save impulser response

%% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot(1:N,h);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')

% Magnitude response
subplot(2,1,2)
dfth = fft(h);
dfth = dfth(1:N/2+1);
plot(0:2*fs/N:fs,mag2db(abs(dfth))); %anders omdat hier fft wordt gebruikt?
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')


%% Filtered white noise Vs. recorded white noise.
sig = wgn(N,1,0); % Generate white noise signal    

% Convolve the IR with the white noise
gen_noise = fftfilt(h,sig);                  

% Send the white noise across the channel
[simin,nbsecs,fs] = initparams(sig,fs,3);
sim('recplay');
rec_noise=simout.signals.values(:,1);

%% Spectrogram and PSD plot
% Spectrogram
figure; 
title('spectogram')
subplot(2,1,1)
spectrogram(sig,N); % Generated white noise
subplot(2,1,2)
spectrogram(rec_noise,N); % Recorded white noise

% PSD
[s1,f1,t1] = spectrogram(sig,N); % Generated white noise
[s2,f2,t2] = spectrogram(rec_noise,N); % Recorded white noise

PSD_gen_noise = sum(abs(s1).^2,2)*fs; % Generated white noise PSD
PSD_rec_noise = sum(abs(s2).^2,2)*fs; % Recorded white noise PSD (ideally only when the signal is active)

figure;
subplot(2,1,1)
plot(f1,pow2db(PSD_gen_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Generated white noise.')
% Output signal
subplot(2,1,2)
plot(f2,pow2db(PSD_rec_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Recorded white noise')
