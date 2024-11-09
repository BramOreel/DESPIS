% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
n = bitsPerPixel; % hier 8
M = 2^n; %2^8 = 256

% QAM modulation
qamStream = qam_mod(bitStream, M); %wat eruit komt is iets van: length(bitStream)/n (idg 19200)

% OFDM modulation
N = length(qamStream); %aantal substreams
Lcp = 16;
ofdmStream = ofdm_mod(qamStream, N, Lcp, 4 );

% Channel
SNR = 10^10;
H = [];
a = 2;
for i = 1:N
    H(i) = 1/(a*i); %channel order = L = 2
end

h = ifft(H');

%{
figure
stem(1:length(h),h)
title('h')
%}


rxOfdmStream = conv(h,ofdmStream); %signaal wordt vervormd door kanaal

% OFDM demodulation
% met channel equalizer
%function [ data_seq, CHANNELS ] = ofdm_demod(OFDM_seq,N,Lcp,varargin, streamLength,channel,equalization )
rxQamStream = ofdm_demod(rxOfdmStream,N,Lcp,M,h,0,1);

% QAM demodulation
rxBitStream = qam_demod(rxQamStream, bitsPerPixel, length(bitStream),4);

% Compute BER
% berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
