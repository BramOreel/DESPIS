%a file containing the experiments needed to run for the first milestone

%% Cleanup
clear all; close all;

% Initialize script parameters
fs = 16000; % Sampling frequency [Hz]
duration = 1;
N = fs*duration; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) used to plot the spectrogra [samples]
windowS = 8000;
Noverlap = windowS/2; 
delay = 200;


%Sinewave
t = linspace(0,2,N);
f_0 = 400;
sinewave = 5*sin(2*pi*f_0*t)';




white_n = wgn(N,1,0); %power is approximately 1 watt, which is 0 dBW

sig = white_n;

% Play and record.
% Call to initparams()
[ simin,nbsecs,fs] = initparams(sig,fs);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);

%Trim the output so we only have the sine
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

yOnset = x; % Determine start of recorded signal [samples]

y = out(yOnset-delay:yOnset + N-delay-1); % Extract the relevant output signal

%Spectrogram of the input and output signal (Welchs method)
%spectrogram(xc,win,noverlap,FFT_LENGTH,fs,'yaxis')

%spectrogram output Welch
[s1,f1,t1,P1] = spectrogram(y,windowS,Noverlap,N/2,fs);  %N/2 moet mss gwn windowS zijn
PSD_rec_noise = sum(P1,2)*fs;

%spectrogram input Welch
[s2,f2,t2,P2] = spectrogram(sig,windowS,Noverlap,N/2,fs);
PSD_gen_noise = sum(P2,2)*fs;

%PSD Bartlett input
[s3,f3,t3,P3] = spectrogram(sig,windowS,0,N/2,fs);
PSD_Bartlett_input = sum(P3,2)*fs;

%PSD Bartlett output
[s4,f4,t4,P4] = spectrogram(y,windowS,0,N/2,fs);
PSD_Bartlett_output = sum(P4,2)*fs;


% Compute and plot the spectrogram

% Input signal
figure; subplot(2,2,1)
%spectrogram(xc,win,noverlap,FFT_LENGTH,fs,'yaxis')
spectrogram(sig,windowS,Noverlap,N/2,fs);
title('Input signal Welch.')

% Output signal
subplot(2,2,2)
spectrogram(y,windowS,Noverlap,N/2,fs);
title('Output signal Welch.')


%spectrogram(xc,win,noverlap,FFT_LENGTH,fs,'yaxis')
subplot(2,2,3)
spectrogram(sig,windowS,0,N/2,fs);
title('Input signal Bartlett.')

% Output signal
subplot(2,2,4)
spectrogram(y,windowS,0,N/2,fs);
title('Output signal Bartlett.')
sgtitle('Spectrogram estimate')





%PSD plots for Welch

figure; subplot(2,1,1)
plot(f2,pow2db(PSD_gen_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')

% Output signal
subplot(2,1,2)
plot(f1,pow2db(PSD_rec_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Welch''s method')

%PSD plots for Bartlett
figure; subplot(2,1,1)
plot(f3,pow2db(PSD_Bartlett_input));
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')
% Output signal
subplot(2,1,2)
plot(f4,pow2db(PSD_Bartlett_output));
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Bartlett''s method')





%% %Session 2; IR1
 


diracsig = zeros(N,1);   %Dirac input signal
diracsig(1) = 1;
sig2 = diracsig;


%Play and record.
% Call to initparams()

[simin,nbsecs,fs] = initparams(sig2,fs,3);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out2=simout.signals.values(:,1);

%Calculate the IR by trimming the output

c2 = max(out2); %impuls respons vinden
x2 = 0;
tol = 10^-9;

for i = 1:size(out2,1)
    if(abs(out2(i) -c2) < tol)
        x2 = i;
        break
    end
end
h = out2(x2-25:x2+150);
% Time domain signal
figure; subplot(2,1,1)
plot(h);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
% Magnitude response
subplot(2,1,2)
plot((0:fs/length(h):fs-fs/length(h)), mag2db(abs(h)));
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')


%% IR2
uMatrix = toeplitz(sig'); % Toeplitz matrix van de ruis
h2 =  lsqr(uMatrix,y,0.1)       ; % Estimate impulse response
%Trim the impulse response
c3 = max(h2); %impuls respons vinden
x3 = 0;

for i = 1:size(h2,1)
    if(abs(h2(i) -c3) < tol)
        x3 = i;
        break
    end
end
h2 = h2(x3-25:x3+150);

save('channel.mat','h2'); % Save impulser response

%% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot(1:N,h2);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')

% Magnitude response
subplot(2,1,2)
dfth = fft(h2);
dfth = dfth(1:N/2+1);
plot(0:2*fs/N:fs,mag2db(abs(dfth))); %anders omdat hier fft wordt gebruikt?
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')


%% Advanced shannon
%run the compute_shannon file

