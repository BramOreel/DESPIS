% Training frames for channel estimation and equalization 
%% Cleanup
clear; clc; close all;

%% Parameters
fs = 16000; % Sampling frequency [Hz]
N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300; % Cyclic prefix length [samples]
M = 16; % QAM constellation size
SNR = 40; % SNR of transmission [dB]
accoustic_transmission = ; % If 1 acoustic transmission occurs, if 0 a simulated transmission.

%% Construct train block.
Nq = log2(M);
train_bits = random([0 1],Nq*(N/2-1),1 ); % Generate a random vector of bits corresponding to a single OFDM frame
%trainblock: (N/2-1)x1
train_block = qam_mod(train_bits,M); % QAM modulate that frame
train_stream = repmat(train_block,100,1); % Repeat the training QAM frame
Tx = ofdm_mod(train_stream,N,Lcp); % OFDM modulate the QAM stream

%% Transmit train block.
if ~accoustic_transmission % Exercise 5.1
    h = load('channel_session5.mat');
    aligned_Rx = fftfilt(h,Tx);
    aligned_Rx = awgn(aligned_Rx,SNR,"measured");
else % Exercise 5.2
    [] = initparams();
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(); % Align input and output
end

%% OFDM Demodulate
[] = ofdm_demod();

%% QAM Demodulate
rx_bits = qam_demod();

%% BER
BER = ber()

%% Plot (real and) estimated channel.
if ~accoustic_transmission
    figure;
    subplot(2,1,1)
    plot();
    title('Real impulse response')
    xlabel('')
    ylabel('')
    subplot(2,1,2)
    plot();
    title('Real frequency response')
    xlabel('')
    ylabel('')    
end

figure;
subplot(2,1,1)
plot();
title('Estimated impulse response')
xlabel('')
ylabel('')
subplot(2,1,2)
plot();
title('Estimated frequency response')
xlabel('')
ylabel('')