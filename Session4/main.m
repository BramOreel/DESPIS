% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');


Nq = bitsPerPixel; % hier 8
M = 2^Nq; % QAM constellation size
N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300;
SNR = 1000000000000000000;

% QAM modulation
[qamStream,x] = qam_mod(bitStream, M);


scatterplot(qamStream);
title('after QAM-modulation')


% OFDM modulation
ofdmStream = ofdm_mod(qamStream, N, Lcp, 4 );

% Channel
L_channel = length(ofdmStream);

%% Channel option 1
%{
h = [];

%impuls respons is van lengte ofdmStream

for i = 1:7
    h = [h; i];
end

%}


h = [1; 0.25; 0.3; 0.9];

%{

min_val = 0;
max_val = 1;
h = [1 randn(1,L_channel-1)];
%h = load('channel_session4.mat').h';
%h = [h];
%h = [1 5 4 3 2 4];

%}

rxOfdmStream = conv(h,ofdmStream);
%rxOfdmStream = awgn(rxOfdmStream,SNR);

% OFDM demodulation
% met channel equalizer

%function [ data_seq, CHANNELS ] = ofdm_demod(OFDM_seq,N,Lcp,varargin, streamLength,channel,equalization )
[ rxQamStream, CHANNELS ] = ofdm_demod(rxOfdmStream,N,Lcp,length(qamStream),h,0,0);

title('after OFDM-modulation')
scatterplot(rxQamStream);

% QAM demodulation
rxBitStream = qam_demod(rxQamStream,M, length(bitStream),x);

% Compute BER
berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
