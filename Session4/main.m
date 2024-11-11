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

figure
plot(1:length(h),h)
title('h')

fs = 16000;
Ntot = length(h);

padLength = abs(mod(size(h,1),N) - N) ;
if(padLength == N)
    padLength = 0;
end

h = [h; zeros(padLength,1)];

s1 = fft(h, N); % H
PSD_h = sum(abs(s1).^2,2)*fs/N;
y = pow2db(PSD_h);

figure;
subplot(2,1,1)
plot(1:N,y);
xlabel('Frequency (kHz)');
ylabel('Normalized Power/frequency (dB/Hz)')
title('|H|^2')

p = length(ofdmStream)/(N+Lcp); %aantal kanalen
fs = 16000; % Sampling frequentie van h
PSD_threshold = 0.5*max(y);


%masker opstellen
f1 = 1:N;
n_off=[]; %alle samples waar masker 0 moet zijn
c = N/fs; % fs/N = one frequency bin

%samples zoeken waarbij vermogen te laag is
for i = 1: length(f1)
    if y(i) <= PSD_threshold
        ni = floor(f1(i)*c);
        a = length(n_off);
        if a ~= 0 
            %geen beginwaarde
            if n_off(a) == ni | ni == 0 %geen DC sample
                n_off = n_off;
            else
                n_off = [n_off ni]; %alle samples in f-dom met te lage vermogen
            end
        else
            %beginwaarde
            if ni ~= 0
                n_off = [n_off ni];
            end
        end
    end
end

MASK = ones(N/2-1,1);
for i = 1:length(n_off)
    MASK(n_off(i),1) = 0;
end

%{
% Gemiddelde vermogen signaal per kanaal
rxOFDM_matrix = reshape(rxOfdmStream,p,[]);
S = [];
for i = 1:p
    S(i) = sum(abs(rxOFDM_matrix(i)).^2)/length(rxOFDM_matrix(i));
end
%}

%{
V=[]; %lengte px1
%Selecteer frequenties met voldoende hoge vermogen = 1; frequenties met slechte vermogen = 0; 
for i=1:length(S)
    if S(i) >= PSD_threshold
        V = [V 1];
    else
        V = [V 0];
    end
end
%}

%{

% Data simuleren en ON-OFF Bit Loading toepassen

data = toeplitz([V(1) fliplr(V(2:end))], V);
data = data(:,1)';
padlength = abs(mod(size(data,2),N/2-1)-(N/2-1));
if padlength == N/2-1
    padlength = 0;
end
data = [zeros(1,size(data,2)); data ; zeros(padlength,size(data,2))];


[qamStream2,x2] = qam_mod(bitStream, M);

ofdmStream2 = ofdm_mod( qamStream2, N, Lcp);

rxOfdmStream2 = fftfilt(h,ofdmStream2);
rxOfdmStream2 = awgn(rxOfdmStream2,SNR);

[ rxQamStream2, CHANNELS2 ] = ofdm_demod(rxOfdmStream2,N,Lcp,length(qamStream2),h,data(:,1),1);

%}

%%

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
scatterplot(qamStream)
title('after QAM_mod')
