% Training frames for channel estimation and equalization 
%% Cleanup
clear; clc; close all;

h = load('channel_session5.mat').h;

%% Parameters

N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300; % Cyclic prefix length [samples]
M = 64; % QAM constellation size
SNR = inf; % SNR of transmission [dB]

accoustic_transmission = 0; % If 1 acoustic transmission occurs, if 0 a simulated transmission.
pilot_tones = 1;            % If 1 pilot tones are used, if 0 the full channel is estimated





%% Construct train block.
Nq = log2(M);
train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits corresponding to a single OFDM frame
%trainblock: (N/2-1)x1
train_block = qam_mod(train_bits,M); % QAM modulate that frame
train_stream = repmat(train_block,10,1); % Repeat the training QAM frame
if pilot_tones

    [Tx,nbOFDMsymb] = ofdm_mod_pilots(train_stream,N,Lcp,train_block);


else

    Tx = ofdm_mod(train_stream,N,Lcp); % OFDM modulate the QAM stream


end
scatterplot(train_stream);
streamlength = length(train_stream);

%% Transmit train block.
if ~accoustic_transmission % Exercise 5.1
    
    aligned_Rx = fftfilt(h,Tx);
    aligned_Rx = awgn(aligned_Rx,SNR,"measured");
else % Exercise 5.2   
    % zelf gekozen hoe lang de synchronisatiepuls is
    Fs = 16000; % Sampling frequency [Hz]
    T = 1/Fs;
    f0 =  100;
    f1 = 2000;
    duration = 1; %Duration of the pulse

    t = 0:T:duration- T;
    pulse = chirp(t,f0,duration,f1)'; %1kHz sin
    sync_pulse = [pulse; zeros(Fs,1)];
    

    %The toplay signal is a sine at another frequency

    toplay = Tx;
    [simin,nbsecs,fs] = initparams(toplay, Fs, sync_pulse);
    sim('recplay');
    Rx = simout.signals.values(:,1); 
    aligned_Rx = alignIO(Rx, fs); % Align input and output

    %This section was just done to be able to receive our training data
    %block. This block has gone through the channel (the air) and has been
    %perturbed by noise. Using this information we can estimate our channel
    %at the output. No actual data has been transmitted here



end

%% OFDM Demodulate
if pilot_tones

    [qam_seq,CHANNELS] = ofdm_demod_pilots(aligned_Rx,N,Lcp,streamlength,train_block,nbOFDMsymb);

else
    [qam_seq,CHANNELS] = ofdm_demod(aligned_Rx,N,Lcp,streamlength,ones(1,N/2-1),train_block);
end





scatterplot(qam_seq);
%Mirror the channels block



%% QAM Demodulate
rx_bits = qam_demod(qam_seq,M,length(train_bits));

%% BER
BER = ber(rx_bits,train_bits);


CHANNELS = [0;CHANNELS ;0; flipud(conj(CHANNELS))];
h_est = ifft(CHANNELS,N);

%% Plot (real and) estimated channel.

    figure;
    subplot(2,1,1)
    plot(h);
    title('Real impulse response')
    xlabel('')
    ylabel('')
    subplot(2,1,2)
    plot(abs(fft(h,N)).^2);
    title('Real frequency response')
    xlabel('')
    ylabel('')    

figure;
subplot(2,1,1)
plot(h_est);
title('Estimated impulse response')
xlabel('')
ylabel('')
subplot(2,1,2)
plot(abs(CHANNELS).^2);
title('Estimated frequency response')
xlabel('')
ylabel('')