%% Cleanup
clear; clc; close all;

%% Parameters
L = 10000; % Scalar length of binary sequence [samples]
SNRs = linspace(1,25,100); % List of SNR values to consider [dB]
N = 6; % List of number of bits per QAM symbol to consider

%% Calculate BER for each constellation-SNR combination.
BERs = [] ; % Placeholder for BERs
for n= 1:N% Loop across values of N
    BERsN = [] ; % Placeholder for BERs
    for SNRidx = 1:size(SNRs,2)% Loop across SNR values
        % Initialize simulation parameters.
        M = 2^n; % QAM constellation size
        SNR = SNRs(SNRidx); % SNR
        
        % Generate a pseudo random binary sequence of a user defined length.
        bit_seq = randi([0, 1], L,1) ;
        
        % Modulate bit sequence.
        QAM_seq = qam_mod(bit_seq,M) ;
        
        % Add white Gaussian noise.
        req_QAM_seq  = awgn(QAM_seq,SNR,"measured");        

        % Demodulate QAM sequence.
        rec_bit_seq = qam_demod(req_QAM_seq,M,L);
        
        % Calculate BER.
        BERsN = [BERsN ;ber(bit_seq,rec_bit_seq)];
    end
    BERs = [BERs BERsN];
end

% Plot results.
[X, Y] = meshgrid([1 2 3 4 5 6],SNRs); % Create a grid of x and y values
surf(X, Y, BERs); % Plot surface
xlabel('Constellation size');
ylabel('SNR');
zlabel('BER');
title('3D Surface Plot of BER vs. SNR and Constellation size');


%semilogy();
%ylabel();
%xlabel();