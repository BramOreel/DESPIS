% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;


% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

% QAM modulation
qamStream = qam_mod(bitStream, bitsPerPixel); 

% OFDM modulation
ofdmStream = ofdm_mod(qamStream, );

% Channel
rxOfdmStream = ofdmStream;

% OFDM demodulation
rxQamStream = rxOfdmStream;

% QAM demodulation
rxBitStream = rxQamStream;

% Compute BER
% berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
