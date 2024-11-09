%% Cleanup
clear; close all; clc;

%% Parameters

Nq = 4; %length of bit sequence
M = 2^Nq; % QAM constellation size
%channel = ; % Impulse response of channel
N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300;
SNR = 1000;

%% Channel effect experiment
% bit stream 

bitStream = randi([0, 1], 1,Nq*(N/2-1))'; % Generate the bit sequence corresponding to 1 QAM symbol

%bitStream = repmat(bitStream,N/2-1,1); % Repeat this bit sequence to fill 1 OFDM frame

% QAM modulation
[qamStream,x] = qam_mod(bitStream,M); % Should be of length N/2-1X1. Each symbol should be the same

% QAM constellation visualization
scatterplot(qamStream); % 

% OFDM modulation
ofdmStream = ofdm_mod(qamStream,N,Lcp,4); % Should be of length N+LcpX1

% Channel
h0 = 1 ;
h1= 0.2;
h2= 0.5;
h = [h0 h1 h2];
rxOfdmStream = fftfilt(h,ofdmStream);
rxOfdmStream = awgn(rxOfdmStream,SNR);

%We should maybe normalise the power again coming from this stream



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


%Why does the cyclic prefic need to be longer than the length of the IR?
%The CP length should be greater than the length of the channel impulse response to ensure that each OFDM symbol can be treated independently without interference from the previous symbol. This allows the channel convolution to be converted into a simple multiplication in the frequency domain.



%Question 5. the filter basically multiplies the whole signal with one in the
%frequency domain. This is a convolution in the time domain with a dirac at
%zero. This shouldnt have any effect. therefore the constellations are
%exactly the same.

%Question6: the constellation scales with the same factor

%Question 7: changing the filter to have a delay strangely ruins
%everything. You get repetions in a circular pattern in the receiver. N/2-1
%repetitions to be exact. The delay in time results in a multplication with
%e^-j2pifktk. this is a phase shift which scales lineairly with f. Higher f
%means more shift. The QAM constellation points are therefore rotated
%around the origin.

%The constant term of one makes it so the constellation points rotate
%around there supposed QAM point. One value for a delay results in a
%linear phase shift. more delay scaling factors result in an amplitutde
%distortion and a spiral pattern

%Question 8:
%Where adding noise on the ideal channel had nearly no effect, now the
%noise highly distorts the spiral pattern and the constellation.

%Question 9: Compensating the channel is highly important!



