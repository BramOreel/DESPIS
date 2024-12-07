% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

% Parameters.
Lh = 200; % Length of impulse response
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
M = 16; % QAM constellation size.
Nq = log2(M);
Lcp = 300; % Cyclic prefix length [samples].
Lt = 10; % Number of training frames
Ld = 0;
Equalization = "adaptive"; % Equalization mode (see ofdm_demod_stereo.m)
mu = 0.02;% NLMS stepsize
alpha = 1; % NLMS regularization
SNR = inf; % SNR of transmission 

% Generate two random impulse responses, and calculate frequency response.
h1 = rand(1,Lh); h2 = rand(1,Lh); % Impulse responses
H =[fft(h1,N/2-1); fft(h2,N/2-1)]; % N/2-1X2 matrix containing frequency transform of h1 and h2

[a,b] = fixed_transmitter_side_beamformer(h1,h2,N); % to do this, we actually need an estimation of our channel first


%% Construct QAM symbol stream.
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M); %Our qam sequence we will use during the adaptive filter
streamLength = length(bitStream);

% Construct train block.
train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits
train_frame = qam_mod(train_bits,M); % QAM training symbols

%We will eventually have a channel estimation here by transmitting training
%frames first

%Now we need a channel estimation for modulation and demodulation

%% OFDM modulation.
[Tx,a,b,nbPackets] = ofdm_mod_stereo(qamStream,N,Lcp,H,Lt,Ld,train_frame,"fixed");

%% Transmit symbol.
Rx = conv(Tx(:,1),h1) + conv(Tx(:,2),h2);
Rx = awgn(Rx,SNR);

%% OFDM demodulation.
%this is just one signal now
[rec_qamStream, CHANNELS] = ofdm_demod_stereo(Rx,N,Lcp,train_frame,Lt,Ld,M,nbPackets,"fixed",mu,alpha);

%% QAM demodulation.
rx_bits = qam_demod();

%% Calculate BER
BER = ber( )
%}