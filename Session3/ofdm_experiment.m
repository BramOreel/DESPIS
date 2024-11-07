%% Cleanup
clear; clc; close all;

%% Parameters
SNR = 30; % Signal-to-noise ratio [dB]
M = 16; % QAM constellation size
N = 16; % Number of symbols per OFDM frame, i.e., the DFT size
L = 10000; % Binary sequence length [samples]
Lcp = 3; % Cyclic prefix length [Samples] (you can ignore this until exercise 3.2.4)

%% OFDM experiment
% Generate a pseudo random binary sequence of a user defined length.
bit_seq = randi([0, 1], 1,L)';

% Modulate bit sequence.
[QAM_seq,x] = qam_mod(bit_seq,M);

% Modulate QAM sequence using OFDM.
OFDM_seq = ofdm_mod(QAM_seq,N,Lcp,4);

% Add white Gaussian noise.
rec_OFDM_seq  = awgn(OFDM_seq,SNR); 

% Demodulate OFDM sequence.
rec_QAM_seq = ofdm_demod(rec_OFDM_seq,N,Lcp,4);

% Demodulate QAM sequence.
rec_bit_seq = qam_demod(rec_QAM_seq,M,L,x);

% Calculate BER.
BER = ber(bit_seq,rec_bit_seq);




%OFDM: first you modulate your signal using QAM, signal is in frequency
%domain. The conversion to ofdm basically means that you put one QAM on one
%OFDM carrier wave. these carriers are at distinct frequencies within the
%bandwith
% 
% 
% 
% Using STFT we can convert a small part of these QAM symbols by 
%doing a ifft where DC = 0 and fn = 0. Imagine the STFT graph.
%The ifft results in a time domain graph of overlapping ofdms?
% 
% the signals overlap but at the right sampling times there is no overlap
% 
% DC = zero because this is a constant power consumption, the bias of the
% signal is also stored hereù
%
%N/2 = 0: we want to avoid spectral leakage
%
%Symmetry because we want ifft to be real
















