%% Cleanup
clear; clc; close all;

%% Parameters
L = 1024; % Scalar length of binary sequence [samples]
SNRs = ; % List of SNR values to consider [dB]
N = 6; % List of number of bits per QAM symbol to consider

%% Calculate BER for each constellation-SNR combination.
BERs = 0 ; % Placeholder for BERs
for n= 1:N% Loop across values of N
    for SNRidx = 1% Loop across SNR values
        % Initialize simulation parameters.
        M = 2^N; % QAM constellation size
        SNR = 10^5; % SNR
        
        % Generate a pseudo random binary sequence of a user defined length.
        bit_seq = randi([0, 1], 1,L)' ;
        
        % Modulate bit sequence.
        QAM_seq = qam_mod(bit_seq,M) ;
        
        % Add white Gaussian noise.
        req_QAM_seq  = awgn(QAM_seq,SNR);        

        % Demodulate QAM sequence.
        %rec_bit_seq = ;
        
        % Calculate BER.
        %BERs() = ;
    end
end

% Plot results.
scatterplot(req_QAM_seq)

%semilogy();
%ylabel();
%xlabel();