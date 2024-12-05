% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters
% Exercise session 7: DMT-OFDM transmission with DD LMS channel tracking.
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300; % Cyclic prefix length [samples].
M = 16; % QAM constellation size.
Nq = log2(M);
SNR = 30; % SNR of transmission [dB].
Lt = 10; % Number of training frames.
fs = 16000; % Sampling frequency [Hz].
channel = "simulation"; % simulation or acoustic

mu = 0.02; % NLMS stepsize
alpha = 1; % NLMS regularization
type = "nlms";
Nswitch = (Lt)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

% bitloading
ON_OFF_mask = ones(N/2-1,1); % Default all bins to one for regular transmission
BW_usage = 60; % Fraction of bins to use for on-off bitloading
bitloading_flag = 0; % If 1 then on-off bitloading is enabled.

%% Determine bit loading
% Channel estimation.
% Construct train block.
train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits
train_block = qam_mod(train_bits,M); % QAM training symbols

%skip this section for now
%{ 

if bitloading_flag
    % sending only trainblock through the channel to estimate channel    
    trainStream = ofdm_mod( );
    
    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session7.mat',smoothing_factor); 
        aligned_Rx = awgn(aligned_Rx,SNR,'measured');
    else
        [] = initparams();
        size(simin)
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO();
    end
    
    % Compute the channel power    
    [~, CHANNELS] = ofdm_demod( );

    % Only keep bins of high energy
    ON_OFF_mask() = ;
end

%}

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M); %Our qam sequence we will use during the adaptive filter
streamLength = length(bitStream);

%% OFDM modulation.
Tx = ofdm_mod(qamStream,N,Lcp,Lt,train_block,ON_OFF_mask);

%% Transmit OFDM sequence.
if channel == "simulation"
    aligned_Rx = simulate_channel(Tx, Nswitch,'channel_session7.mat',smoothing_factor);
    aligned_Rx = awgn(aligned_Rx,SNR,'measured');
elseif channel == "acoustic"
    duration = 1;
    T = 1/fs;
    t = 0:T:duration- T;
    pulse = chirp(t,100,duration,2000)';
    sync_pulse = [pulse; zeros(fs,1)];

    [simin,nbsecs] = initparams(Tx,fs,sync_pulse);
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(Rx,fs);
end

%% OFDM Demodulate
[rx_qam, CHANNELS] = ofdm_demod(aligned_Rx,N,Lcp,Lt,M,train_block,ON_OFF_mask,mu,alpha,type);

%% QAM demodulation.
rx_bits = qam_demod(rx_qam,M,streamLength);

%% Bit error rate
BER = ber(rx_bits,bitStream);

% Construct image from bitstream
imageRx = bitstreamtoimage(rx_bits, imageSize, bitsPerPixel);

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
disp(BER);