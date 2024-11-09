% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

% QAM modulation
qamStream = qam_mod(bitStream, bitsPerPixel); 

% OFDM modulation
ofdmStream = ofdm_mod(qamStream, 16, 7,4 );

% Channel
SNR = 10^10;
H = [];
a = 2;
for i = 1:length(ofdmStream)
    H(i) = 1/(a*i); %channel order = L = 2
end
h = ifft(H');

rxOfdmStream = conv(h,ofdmStream);

% OFDM demodulation
rxQamStream = ofdm_demod(rxOfdmStream, 16, 7,4);

% QAM demodulation
rxBitStream = qam_demod(rxQamStream, bitsPerPixel, length(bitStream),4);

% Compute BER
% berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);

% Plot images
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
