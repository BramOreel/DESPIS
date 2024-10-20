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
y(1) = 1;




sig = y;                       

%% Play and record.
% Call to initparams()


[simin,nbsecs,fs] = initparams(sig,fs,3);
% Call to recplay.mdl to play simin and record simout
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);



%% Calculate the IR by trimming the output

c = max(out); %impuls respons vinden
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
plot(1:N, h);
xlabel('Time [samples]')
ylabel('Impulse response [arb.]')
% Magnitude response
subplot(2,1,2)
plot((0:fs/length(h):fs-fs/length(h)), mag2db(abs(h)));
xlabel('Frequency [Hz]')
ylabel('Magnitude response [dB]')




%% 
%{
First, use IR1.m to make a new estimate of the IR. 


Then,
without moving the microphone, repeat the white-noise experiment from
exercise 1-2. 


If you now convolve the white noise signal with the estimated
IR (use the command fftfilt), this should yield an output signal with
similar characteristics as the recorded signal (why?). 

%}


% Filtered white noise Vs. recorded white noise.
sig = wgn(N,1,0); % Generate white noise signal
[simin,nbsecs,fs] = initparams(sig,fs,3);
sim('recplay');
unfiltered_noise= simout.signals.values(:,1);


%filtered noise

i = 1;
max_signal = max(abs(unfiltered_noise))

filtered_noise = unfiltered_noise;

while i <= length(filtered_noise)
    if(abs(filtered_noise(i)) >= 0.2*max_signal) %deze waarde moet eruit
            filtered_noise(i) = 0;
            i = i - 1;
    end
    i = i + 1;
end


figure; subplot(4,1,1)
plot(unfiltered_noise)
ylim([-0.08 0.08])
ylabel('unfiltered noise')

subplot(4,1,2)
plot(filtered_noise)
ylim([-0.08 0.08])
ylabel('filtered noise')




% Convolve the IR with the white noise
tot_noise = fftfilt(h,filtered_noise); %IR convolueren met gefilterde ruis signaal
%this should yield an output signal with similar characteristics as the recorded signal (why?).

%dirac in frequentiedom is ook een dirac

% Send the white noise across the channel
[simin,nbsecs,fs] = initparams(tot_noise,fs,3);
sim('recplay');
rec_noise= simout.signals.values(:,1);

subplot(4,1,3)
plot(tot_noise)
ylabel('h*noise')

subplot(4,1,4)
plot(rec_noise)
ylabel('recorded convolution: h*noise') %het is in tijd verschoven


title('Noise')


%% Spectrogram and PSD plot
% Spectrogram

%{
Compare the spectrograms and PSDs of the recorded signal and the convolved signal. Do
they look more or less the same? For the PSD estimation, you can use
Welch’s method.
%}


%spectrograms
figure; 
title('spectogram')
subplot(2,1,1)
spectrogram(tot_noise,N); % Generated white noise
subplot(2,1,2)
spectrogram(rec_noise,N); % Recorded white noise

% PSD
[s1,f1,t1] = spectrogram(sig,N); % Generated white noise
[s2,f2,t2] = spectrogram(rec_noise,N); % Recorded white noise

PSD_gen_noise = sum(abs(s1).^2,2)/T; % Generated white noise PSD
PSD_rec_noise = sum(abs(s2).^2,2)/T; % Recorded white noise PSD (ideally only when the signal is active)

figure;
title('PSD')
subplot(2,1,1)
plot(f1,pow2db(PSD_gen_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
subtitle('Generated white noise.')

% Output signal
subplot(2,1,2)
plot(f2,pow2db(PSD_rec_noise));
xlabel('Frequency (kHz)');
ylabel('Power/frequency (dB/Hz)')
subtitle('Recorded white noise')




