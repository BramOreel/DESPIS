% Training frames for channel estimation and equalization 
%% Cleanup
clear; clc; close all;

%% Parameters
fs = 16000; % Sampling frequency [Hz]
N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300; % Cyclic prefix length [samples]
M = 16; % QAM constellation size
SNR = 40; % SNR of transmission [dB]

accoustic_transmission = 1; % If 1 acoustic transmission occurs, if 0 a simulated transmission.

%% Construct train block.
Nq = log2(M);
train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits corresponding to a single OFDM frame
%trainblock: (N/2-1)x1
train_block = qam_mod(train_bits,M); % QAM modulate that frame
train_stream = repmat(train_block,100,1); % Repeat the training QAM frame
Tx = ofdm_mod(train_stream,N,Lcp); % OFDM modulate the QAM stream
streamlength = length(train_stream);

%% Transmit train block.
if ~accoustic_transmission % Exercise 5.1
    h = load('channel_session5.mat').h;
    aligned_Rx = fftfilt(h,Tx);
    aligned_Rx = awgn(aligned_Rx,SNR,"measured");
else % Exercise 5.2    
    h = load('channel_session5.mat').h;
    % zelf gekozen hoe lang de synchronisatiepuls is
    pulse = sin(2*pi*400*(0:1/fs:1-1/fs))';
    
    sync_pulse = [pulse; zeros(length(h),1)];
    [simin,nbsecs,fs] = initparams(Tx, fs, sync_pulse);
    sim('recplay');
    Rx = simout.signals.values(:,1); 
    aligned_Rx = alignIO(Rx, fs); % Align input and output
end

%% OFDM Demodulate
[qam_seq,CHANNELS] = ofdm_demod(aligned_Rx,N,Lcp,streamlength,ones(1,N/2-1),train_block);

%Mirror the channels block
CHANNELS = [0;CHANNELS ;0; flipud(conj(CHANNELS))];
h_est = ifft(CHANNELS,500);


%% QAM Demodulate
rx_bits = qam_demod(qam_seq,M,length(train_bits));

%% BER
BER = ber(rx_bits,train_bits);

%% Plot (real and) estimated channel.
if ~accoustic_transmission
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
end

figure;
subplot(2,1,1)
plot(ifft(CHANNELS,N));
title('Estimated impulse response')
xlabel('')
ylabel('')
subplot(2,1,2)
plot(abs(CHANNELS).^2);
title('Estimated frequency response')
xlabel('')
ylabel('')