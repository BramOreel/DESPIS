% Plays and records signals, which are subsequently analyzed by looking at
% the spectrograms, and the power spectral densities.

%% Cleanup
clear all; close all;

%% Initialize script parameters
fs = 16000; % Sampling frequency [Hz]
N = 32000; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) used to plot the spectrogra [samples]
Noverlap = N/2; 

%% Construct signals
t = linspace(0,2,N);
f0 = 400;
f1 = 100;
f2 = 200;
f3 = 500;
f4 = 1000;
f5 = 1500;
f6 = 2000;
f7 = 4000;
f8 = 6000;


f_0 = 1500;

sinewave = 5*sin(2*pi*f_0*t)';
sines = sin(2*pi*f1*t)' + sin(2*pi*f2*t)' + sin(2*pi*f3*t)'+sin(2*pi*f4*t)'+sin(2*pi*f5*t)'+sin(2*pi*f6*t)'+sin(2*pi*f7*t)'+sin(2*pi*f8*t)';

white_n = wgn(N,1,0); %power is approximately 1 watt, which is 0 dBW

sig = sinewave; 

%% Play and record.
% Call to initparams()
[ simin,nbsecs,fs] = initparams(sig,fs);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%% Compute and plot the spectrogram

% Input signal
figure; subplot(2,2,1)
%spectrogram(xc,win,noverlap,FFT_LENGTH,fs,'yaxis')
spectrogram(sig,hamming(N),Noverlap,N,fs,"power")
title('Input signal.')

% Output signal
subplot(2,1,2)
spectrogram(out,hamming(N),Noverlap,N,fs,"power")
title('Output signal.')

%% Select input and output signals to compute PSD
% Ideally the output signal should be trimmed such that only periods where
% the signal is active are considered as the PSD assumes stationarity. 
% Else, the PSD will be biased due to the inclusion of the silence.


in = sig;
out = out;

%% Compute and plot the power spectral density (PSD)...



% ...Using Welch's method
% Input signal



%{

Calculate the PSD by averaging multiple spectrum estimates, all with a DFT size equal to N.
s = spectrogram(x,window,noverlap,nfft, fs) returns the Short-Time Fourier Transform (STFT) of the 
input signal x. Each column of s contains an estimate of the short-term, 
time-localized frequency content of x. The magnitude squared of s is known 
as the spectrogram time-frequency representation of x.
nfft sampling points to calculate the discrete Fourier transform
%}


%{
[s,w,t] = spectrogram(___) returns a vector of normalized frequencies, w, 
and a vector of time instants, t, at which the STFT is computed. 
s = spectrogram(x,window,noverlap,nfft, fs) returns the Short-Time Fourier Transform (STFT) of the 
input signal x. Each column of s contains an estimate of the short-term, 
time-localized frequency content of x. The magnitude squared of s is known 
as the spectrogram time-frequency representation of x.

This syntax can include any combination of input arguments from previous
syntaxes.


https://www.youtube.com/watch?v=YK1F0-3VvQI
%}

segmentLenght = 1000;
noverlap = segmentLenght/2;
[PSD_Welch_input,Fin] = pwelch(in,segmentLenght,noverlap,N,fs);
[PSD_Welch_output,Fout] = pwelch(out,segmentLenght,noverlap,N,fs);




%[Sin,Win,Tin] = spectrogram(in,N,Noverlap,N,fs);
%[Fin,PSD_Welch_input] = pwelch(Sin);%is dit juist? Wrm Sin en niet gwn het signaal 'in'?
% Output signal
%[Sout,Wout,Tout] = spectrogram(out,N,Noverlap,N,fs);
%[Fout,PSD_Welch_output] = pwelch(Sout);






% Plot results


% Input signal
figure; subplot(2,1,1)
plot(Fin,pow2db(PSD_Welch_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')

% Output signal
subplot(2,1,2)
plot(Fout,pow2db(PSD_Welch_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Welch''s method')


% ...Using Bartlett's method
% Input signal
%[BSin,BFin,BTin] = spectrogram(in,rectwin(N),0,N,fs);
[PSD_Bartlett_input, BFin] = pwelch(in,rectwin(N),0,N,fs);
% Output signal
%[BSout,BFout,BTout] = spectrogram(out,rectwin(N),0,N,fs);
[PSD_Bartlett_output, BFout] = pwelch(out,rectwin(N),0,N,fs);
% Plot results

% Input signal
figure; subplot(2,2,1)
plot(BFin,pow2db(PSD_Bartlett_input));
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,1,2)
plot(BFout,pow2db(PSD_Bartlett_output));
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Bartlett''s method')


 
% ...using the magnitude squared of the frequency spectrum
% Input signal
%https://www.youtube.com/watch?v=pfjiwxhqd1M&t=110s
%opm: magnitude van de fft zegt niets over de amplitude van het signaal!



Nin = length(in); %3200
indft = fft(in); %3200
indft = abs(indft).^2; %magnitude squared


Power_fft_squared_input = indft; % Compute magnitude absolute value squared of the fft of the input
% Rescaling and computations to calcualte the one sided PSD to be
% consistent with 'periodogram' (You can ignore the following lines)
Power_fft_squared_input = Power_fft_squared_input/(length(in)*fs);
PSD_fft_squared_input = Power_fft_squared_input(1:floor(length(in)/2)+1);
PSD_fft_squared_input(2:ceil(length(in)/2)) = ...
    PSD_fft_squared_input(2:ceil(length(in)/2)) + ...
    flipud(Power_fft_squared_input(floor(length(in)/2)+2:length(in)));


% Ouput signal
Nout = length(out);
outdft = fft(out);


%outdft = outdft(1:Nout-1);
outdft = abs(outdft).^2;


Power_fft_squared_output = outdft; % Compute magnitude absolute value squared of the fft of the output
% Rescaling and computations to calculate the one sided PSD to be
% consistent with 'periodogram' (You can ignore the following lines)
Power_fft_squared_output = Power_fft_squared_output/(length(out)*fs);
PSD_fft_squared_output = Power_fft_squared_output(1:floor(length(out)/2)+1);
PSD_fft_squared_output(2:ceil(length(out)/2)) = ...
    PSD_fft_squared_output(2:ceil(length(out)/2)) + ...
    flipud(Power_fft_squared_output(floor(length(out)/2)+2:length(out)));

% Plot results
%lengte = length(0:fs/length(in):fs/2) = 16001
%lengte2 = length(Power_fft_squared_output) = 80128

% Input signal
%reeel signaal, dus slechts helft van het spectrum nodig
Power_fft_squared_input = Power_fft_squared_input(1:floor(Nin/2+1));

figure; subplot(2,2,1)

plot((0:fs/length(in):fs/2),pow2db(Power_fft_squared_input));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal

Power_fft_squared_output = Power_fft_squared_output(1:floor(Nout/2+1));
subplot(2,1,2)


plot((0:fs/length(out):fs/2),pow2db(Power_fft_squared_output));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using the magnitude squared of the frequency spectrum')
%}