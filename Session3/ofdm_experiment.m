%% Cleanup
clear; clc; close all;

%% Parameters
SNR = ; % Signal-to-noise ratio [dB]
M = ; % QAM constellation size
N = ; % Number of symbols per OFDM frame, i.e., the DFT size
L = ; % Binary sequence length [samples]
Lcp = ; % Cyclic prefix length [Samples] (you can ignore this until exercise 3.2.4)

%% OFDM experiment
% Generate a pseudo random binary sequence of a user defined length.
bit_seq = ;

% Modulate bit sequence.
QAM_seq = ;

% Modulate QAM sequence using OFDM.
OFDM_seq = ;

% Add white Gaussian noise.
rec_OFDM_seq  = ;

% Demodulate OFDM sequence.
rec_QAM_seq = ;

% Demodulate QAM sequence.
rec_bit_seq = ;

% Calculate BER.
BER = ;