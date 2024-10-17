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
dftsize = 100; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) [samples]
Noverlap = 0; 



%% Create the signal to be played

T = 1/fs;
duration = 2;
N = duration*fs; %totaal aantal samples

y = zeros(N,1);
y = [y;1];
y = [y;zeros(N,1)];



sig = y;                       

%% Play and record.
% Call to initparams()


[simin,nbsecs,fs] = initparams(sig,fs,3);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);



%% Calculate the IR by trimming the output

c = max(out);
x = 0;
tol = 10^-9

for i = 1:size(out,1)
    if(abs(out(i) -c) < tol)
        x = i;
        break
    end
    
end




h = out(x-25:x+150);


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



%% Filtered white noise Vs. recorded white noise.
sig = wgn(N,1,0); % Generate white noise signal    

% Convolve the IR with the white noise
gen_noise = fftfilt(h,sig);                  

% Send the white noise across the channel
[simin,nbsecs,fs] = initparams(gen_noise,fs,3);
sim('recplay');
rec_noise= simout.signals.values(:,1);


%% Spectrogram and PSD plot
% Spectrogram
spectrogram(sig,fs); % Generated white noise
spectrogram(rec_noise,fs); % Recorded white noise

% PSD
[s1,f1,t1] = spectrogram(sig,fs); % Generated white noise
[s2,f2,t2] = spectrogram(rec_noise,fs); % Recorded white noise

PSD_gen_noise = sum(abs(s1).^2,2)/T; % Generated white noise PSD
PSD_rec_noise = sum(abs(s2).^2,2)/T; % Recorded white noise PSD (ideally only when the signal is active)

figure;
subplot(1,2,1)
plot(f1,pow2db(PSD_gen_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Generated white noise.')
% Output signal
subplot(1,2,2)
plot(f2,pow2db(PSD_rec_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Recorded white noise')
