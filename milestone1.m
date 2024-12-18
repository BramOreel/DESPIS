%a file containing the experiments needed to run for the first milestone

%% Cleanup
clear all; close all;

%1. The white noise experiment of exercise 1-2 (week 1).

% Initialize script parameters
fs = 16000; % Sampling frequency [Hz]
duration = 1;
N = fs*duration; % Discrete Fourier Transform (DFT) size [Samples] 
% Overlap length between subsequent frames in the
% short-time-Fourier-transform (STFT) used to plot the spectrogra [samples]
windowS = 1024;
Noverlap = windowS/2; 
delay = 200;


%Sinewave
t = linspace(0,duration,N);
f_0 = 400;
sinewave = sin(2*pi*f_0*t)';

white_n = wgn(N,1,0); %power is approximately 1 watt, which is 0 dBW

sig = white_n;

% Play and record.
% Call to initparams()
[ simin,nbsecs,fs] = initparams(sig,fs,3);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);


%%
%2. The estimation of the IR with IR2.m (use the simout from previous
%experiment without playing a new white noise sequence).

%Trim the output so we only have the noise signal
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

% Extract the relevant output signal
if yOnset-delay <0 %anders krijg je negatieve indices
    filterd_output = out(1:N-1);
else
   filterd_output = out(yOnset-delay:yOnset + N-delay-1);
end


%%  Spectogram

%spectrogram(xc,win,noverlap,FFT_LENGTH,fs,'yaxis')
%Compute and plot the spectrogram

[s2,f2,t2,P2] = spectrogram(sig,hamming(windowS),Noverlap,windowS,fs);
[s1,f1,t1,P1] = spectrogram(filterd_output,hamming(windowS),Noverlap,windowS,fs);  %N/2 moet mss gwn windowS zijn
[s3,f3,t3,P3] = spectrogram(sig,windowS,0,windowS,fs);
[s4,f4,t4,P4] = spectrogram(filterd_output,windowS,0,windowS,fs);


figure; 
subplot(2,2,1)
%spectrogram input Welch
spectrogram(sig,hamming(windowS),Noverlap,windowS,fs);
clim([-200 0])
title('Input signal Welch.')
PSD_gen_noise = sum(P2,2)*fs;

subplot(2,2,2)
spectrogram(filterd_output,hamming(windowS),Noverlap,windowS,fs);
clim([-200 0])
%spectrogram output Welch
title('Output signal Welch.')
PSD_rec_noise = sum(P1,2)*fs;

subplot(2,2,3)
spectrogram(sig,windowS,0,windowS,fs);
clim([-200 0])
%PSD Bartlett input
title('Input signal Bartlett.')
PSD_Bartlett_input = sum(P3,2)*fs;

%PSD Bartlett output
subplot(2,2,4)
spectrogram(filterd_output,windowS,0,windowS,fs);
clim([-200 0])
title('Output signal Bartlett.')
PSD_Bartlett_output = sum(P4,2)*fs;
sgtitle('Spectrogram estimate')

%%  PSD plot

%PSD plots for Welch
figure; subplot(2,1,1)
plot(f2,pow2db(PSD_gen_noise));
ylim([0 30]);
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')

% Output signal
subplot(2,1,2)
plot(f1,pow2db(PSD_rec_noise));
ylim([-65 5])
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Welch''s method')

%PSD plots for Bartlett
figure; subplot(2,1,1)
plot(f3,pow2db(PSD_Bartlett_input));
ylim([0 30]);
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)')
title('Input signal.')

% Output signal
subplot(2,1,2)
plot(f4,pow2db(PSD_Bartlett_output));
ylim([-65 5]);
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)')
title('Output signal.')
sgtitle('Power Spectral Density estimate using Bartlett''s method')



%%
uMatrix = toeplitz(sig'); % Toeplitz matrix van de ruis
h2 =  lsqr(uMatrix,filterd_output,0.1)       ; % Estimate impulse response
%Trim the impulse response
c3 = max(h2); %impuls respons vinden
x3 = 0;
tol = 10^-9;


for i = 1:size(h2,1)
    if(abs(h2(i) -c3) < tol)
        x3 = i;
        break
    end
end
if x3-25 <0 %anders krijg je negatieve indices
    range=1:175;
    h2_trim = h2(range);
else
    range = x3-25:x3+150;
    h2_trim = h2(range);
end



save('IRest.mat','h2_trim'); % Save impulse response

% Plot IR.
% Time domain signal
figure; subplot(2,1,1)
plot((range),h2_trim);
%ylim([-0.01 0.01]);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
title('time response')

% Magnitude response
subplot(2,1,2)
dfth = fft(h2);
dfth = dfth(1:N/2+1);
plot(0:2*fs/N:fs,mag2db(abs(dfth))); %anders omdat hier fft wordt gebruikt? % 
ylim([-100 20]);
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')
title('frequency response')
sgtitle('Impulse response calculqted using lsqr method')



%% The estimation of the IR with IR1.m.
 
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

yOnset2 = x2;
if yOnset2-delay <0 %anders krijg je negatieve indices
    delayed_output = out2(1:N-1);
else
   delayed_output = out2(yOnset2-delay:yOnset2 + N-delay-1);
end

h = out2(x2-25:x2+150);
% Time domain signal
figure; subplot(2,1,1)
plot((x2-25:x2+150),h);
%ylim([-0.01 0.01]);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
title('time response')


% Magnitude response
subplot(2,1,2)
dfth = fft(delayed_output);
dfth = dfth(1:N/2+1);
plot(0:2*fs/N:fs,mag2db(abs(dfth)));
ylim([-100 20]);
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')
title('Frequency response')
sgtitle('Impulse response from dirac recording')




%% Advanced shannon
%run the compute_shannon file

%% 2-3
%run the 2-3 file








