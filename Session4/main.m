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
MASK = eye(N/2-1,1);
BWusage = 0.5;

% QAM modulation
[qamStream,x] = qam_mod(bitStream, M);


scatterplot(qamStream)
title('after QAM_mod')


% OFDM modulation
ofdmStream = ofdm_mod(qamStream, N, Lcp, 4 );

% Channel
h = load('channel_session4.mat').h';
h = h';

%%
rxOfdmStream = fftfilt(h,ofdmStream);
rxOfdmStream = awgn(rxOfdmStream,SNR);

% OFDM demodulation

%function [ data_seq, CHANNELS ] = ofdm_demod(OFDM_seq,N,Lcp,varargin, streamLength,channel,MASK,equalization )
[ rxQamStream, CHANNELS ] = ofdm_demod(rxOfdmStream,N,Lcp,length(qamStream),h,MASK,1);

scatterplot(rxQamStream);
title('after OFDM-modulation')

%%

% QAM demodulation
rxBitStream = qam_demod(rxQamStream,M, length(bitStream),x);

% Compute BER
berTransmission = ber(bitStream,rxBitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rxBitStream, imageSize, bitsPerPixel);
%% OEF 4.3


[qamStream2,x2] = qam_mod(bitStream, M);

ofdmStream2 = ofdm_mod( qamStream2, N, Lcp);

rxOfdmStream2 = fftfilt(h,ofdmStream2);
rxOfdmStream2 = awgn(rxOfdmStream2,SNR);


figure
plot(1:length(h),h)
title('h')

fs = 16000;% Sampling frequentie van h


sh = N/2-1; %aantal samples van H
s1 = fft(h, sh); % H, fft pads h already to length sh

PSD_h = sum(abs(s1).^2,2)*fs/sh;
y = PSD_h; %vermogen in dB

figure;
subplot(1,1,1)
plot(1:sh,y);
xlabel('sample');
ylabel('|H|^2')

%p = length(ofdmStream)/(N+Lcp); %aantal kanalen
fs = 16000; 
PSD_threshold = 0.5*max(abs(y));

%masker opstellen
n_on=[]; %alle samples waar masker 1 moet zijn
c = fs/N; % fs/N = one frequency bin
y_below = [];

%samples zoeken waarbij vermogen hoog is
for i = 1:sh
    if y(i) >= PSD_threshold %vermogen die hoger is dan de treshold
        %ni = floor(i*c); %frequentie bin: k*fs/N
        a = length(n_on);
        if a ~= 0 
            %geen beginwaarde
            if n_on(a) == i | i == 0
                n_on = n_on;

            else
                n_on = [n_on [i; y(i)]];
            end
        else
            %beginwaarde
            if i ~= 0
                n_on = [n_on [i; y(i)]];
            end
        end
    else
        y_below= [y_below y(i)];
    end

end

%BW to use
nuse = sort(n_on(1,:),'descend'); %alle frequency bins
nuse = nuse(1:floor(BWusage*length(nuse)));

MASK = zeros(sh,1);
for i = 1:size(nuse,2)
    MASK(nuse(i),1) = 1;
end


%{
% Gemiddelde vermogen signaal per kanaal
rxOFDM_matrix = reshape(rxOfdmStream,p,[]);
S = [];
for i = 1:p
    S(i) = sum(abs(rxOFDM_matrix(i)).^2)/length(rxOFDM_matrix(i));
end
%}


[ rxQamStream2, CHANNELS2 ] = ofdm_demod(rxOfdmStream2,N,Lcp,length(qamStream2),h,MASK,1);

%}

%%

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
scatterplot(qamStream)
title('after QAM_mod')
