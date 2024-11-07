%% Cleanup
clear; close all; clc;

%% Parameters

Nq = 5; %length of bit sequence
M = 2^Nq; % QAM constellation size
%channel = ; % Impulse response of channel
N = 16; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 4;

%% Channel effect experiment
% bit stream 

bitSeq = randi([0, 1], 1,Nq)'; % Generate the bit sequence corresponding to 1 QAM symbol
bitStream = repmat(bitSeq,N/2-1,1); % Repeat this bit sequence to fill 1 OFDM frame

% QAM modulation
[qamStream,x] = qam_mod(bitStream,M); % Should be of length N/2-1X1. Each symbol should be the same

% QAM constellation visualization
scatterplot(qamStream); % 

% OFDM modulation
ofdmStream = ofdm_mod(qamStream,N,Lcp,4); % Should be of length N+LcpX1

% Channel
h0 =1;
h1=0;
h2=0;
h = [h0 h1 h2];
rxOfdmStream = fftfilt(h,ofdmStream);
%rxOfdmStream = awgn();

% OFDM demodulation
rxQamStream = ofdm_demod(rxOfdmStream,N,Lcp,4);

% QAM constellation visualization
scatterplot(rxQamStream);

% QAM demodulation
rxBitStream = qam_demod(rxQamStream,M,size(bitStream,1),x);

% Compute BER
berTransmission = ber(bitStream,rxBitStream);




% The name Orthogonal Frequency Division Multiplexing (OFDM) comes from two key concepts:
% 
% Orthogonal: Refers to the fact that the subcarriers used in OFDM are mathematically orthogonal to each other. This means the signals on each subcarrier do not interfere with one another, allowing efficient use of bandwidth.
% 
% Frequency Division Multiplexing (FDM): This is a method where multiple signals are transmitted over different frequency bands (or subcarriers) simultaneously. OFDM is a specific form of FDM where the subcarriers are spaced such that they remain orthogonal.
% 
% In short, OFDM is a method of dividing a broad frequency band into multiple narrowband subcarriers that are orthogonal to each other, enabling efficient data transmission.
