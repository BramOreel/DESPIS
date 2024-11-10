% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');


Nq = bitsPerPixel; % hier 8
M = 2^Nq; % QAM constellation size
N = 2048; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300;
SNR = 10000000;

% QAM modulation
[qamStream,x] = qam_mod(bitStream, M);


scatterplot(qamStream)
title('after QAM_mod')


% OFDM modulation
ofdmStream = ofdm_mod(qamStream, N, Lcp, 4 );
%display(size(ofdmStream),'ofdstream')

% Channel
h = load('channel_session4.mat').h';
h = [h'];

figure
plot([1:length(h)],h)
title('h')
%{
fs = 16000;
Nh = length(h);
segmentLength = length(h)/30;
noverlap = segmentLength/2;
[PSD_Welch_input,Fin] = pwelch(h,segmentLength,noverlap,1/fs,fs);

%}
%{
for i =1 : length(ofdmStream)
    h(i) = 1/i;
end

%}


rxOfdmStream = fftfilt(h,ofdmStream);
%rxOfdmStream = awgn(rxOfdmStream,SNR);

% OFDM demodulation

%function [ data_seq, CHANNELS ] = ofdm_demod(OFDM_seq,N,Lcp,varargin, streamLength,channel,MASK,equalization )
[ rxQamStream, CHANNELS ] = ofdm_demod(rxOfdmStream,N,Lcp,length(qamStream),h,0,1);

h_hat = ifft(CHANNELS)


scatterplot(rxQamStream);
title('after OFDM-modulation')

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
