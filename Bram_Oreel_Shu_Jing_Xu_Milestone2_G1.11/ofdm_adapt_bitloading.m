% DMT-OFDM transmission scheme
%% Cleanup
clear; close all; clc;

% Convert BMP image to bitstream
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');

N = 1024;
Lcp = 500;
T = 10;
SNR = 10;
channel = load('channel_session4.mat').h';
equalization = 1;
streamLength = length(bitStream);

%Now we don't use on off bit loading, rather we vary M, while using every
%frequency
%Determine M for each frequency
H = fft(channel, N);       % Channel frequency response
H_abs = abs(H).^2;

b_mat = floor(log2(1 + (N/2-1)*db2pow(SNR).*H_abs()/(T*sum(H_abs(1:N/2-1),"all"))));
bit_column = sum(b_mat(1:N/2-1),"all");
M_vary = 2.^b_mat;





%QAM modulation will be a more intensive step
%We have to modulate data that will sit on certain frequencies with certain
%constellation sizes
%Just modulate everything and but a zero if M = 1
j = 1; %Where are we in the bitstream?
qamStream = [];
while j <= streamLength
    for Nq = b_mat(1:N/2-1)
    M = 2^Nq;

    if(M == 1) && j < streamLength
        qamStream = [qamStream ; 0];


    elseif M ~=1 && j + Nq - 1 <= streamLength
         qamStream = [qamStream ; qam_mod(bitStream(j:j+Nq-1),M)]; %vreemde N HIER
        j = j+Nq;
    elseif M ~= 1 && j <= streamLength
        qamStream = [qamStream ; qam_mod(bitStream(j:end),M)];
        j = j+Nq;
    end
    end
end


padLength = abs(mod(size(qamStream,1),N/2-1) -(N/2-1)) ; % Number of bits to append such that it can be divided nicely into the M-ary QAM format
if(padLength == (N/2-1))
    padLength = 0;
end
qamStream_padded = [qamStream ;zeros(padLength,1)];


ofdmStream = ofdm_mod(qamStream_padded,N,Lcp);

%Send the ofdm stream through the channel
rxOfdmStream = fftfilt(channel,ofdmStream);
rxOfdmStream = awgn(rxOfdmStream,SNR,'measured');

%ofdm demodulation
rxQamStream = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), channel, ones(1,N/2-1), equalization);


%convert the data back to a bitstream

j  = 1;
bitMatrix = [];
while j <= length(qamStream)
    for Nq = b_mat(1:N/2-1)
    M = 2^Nq;
    if M ~=1 && j == length(qamStream)
        nb_lastBits = streamLength - length(bitMatrix);
        lastBits_padded = qam_demod(rxQamStream(j), M, Nq);
        bitMatrix = [bitMatrix; lastBits_padded(end - nb_lastBits + 1:end)];
    elseif M ~= 1 && j < length(qamStream)
        bitMatrix = [bitMatrix ; qam_demod(rxQamStream(j),M,Nq)];
    end
    j = j+1;

    end
end

berTransmission = ber(bitStream,bitMatrix);

% Construct image from bitstream
imageRx = bitstreamtoimage(bitMatrix, imageSize, bitsPerPixel);

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;







