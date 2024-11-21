% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters.
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300; % Cyclic prefix length [samples]. Lcp has to be bigger than N/2-1 i think
Nq = 4;
M = 2^Nq; %  constellation size.
SNR = 15; % SNR of transmission [dB].
Lt = 5; % Number of training frames.
Ld = 5; % Number of data frames.
fs = 16000; % Sampling frequency [Hz].
channel = "simulation"; % acoustic or simulation

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);
streamLength = length(bitStream);



% Bitloading
ON_OFF_mask = ones(1,N/2-1); % Default all bins to one for regular transmission
bitloading_flag = 0; % If 1 then on-off/adaptive bitloading is enabled.
bitloading_type = 0; % on-off or adaptive 
BW_usage = 50; % Fraction of bins to use for on-off bitloading
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

%% Determine bit loading
% Channel estimation based on a dummy transmission for bitloading

train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits
train_block = qam_mod(train_bits,M); % QAM modulate

%{
if bitloading_flag

    % Dummy transmission
    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session6.mat',smoothing_factor);
        aligned_Rx = awgn(aligned_Rx,SNR,'measured');
    else
        [x,y] = initparams();
        size(simin)
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO();
    end

    % Extract the estimated channels
    [~, CHANNELS] = ofdm_demod();

    if bitloading_type == "on-off"
        % Only keep bins of high energy
        ON_OFF_mask() = ; % ON-OFF mask with 1 denoting the usage of a bin.
    elseif bitloading_type == "adaptive"
        M = ;     % Constellation sizes
    end
end
%}


%% OFDM modulation
[ Tx, nbPackets ] = ofdm_mod(qamStream,N,Lcp,ON_OFF_mask,Lt,Ld,train_block);

%% Transmit OFDM sequence.
if channel == "simulation"
    aligned_Rx = simulate_channel(Tx, Nswitch,'channel_session6.mat',smoothing_factor);
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


    %Hoe relateert de Lcp aan de channel estimate? N?
end

%% OFDM Demodulate
[rx_qam, CHANNELS] = ofdm_demod(aligned_Rx,N,Lcp,length(qamStream),ON_OFF_mask,train_block,Lt,Ld,nbPackets);

%% QAM Demodulate
rx_bits = qam_demod(rx_qam,M,streamLength);

%% Bit error rate
BER = ber(rx_bits,bitStream );

% Construct image from bitstream
imageRx = bitstreamtoimage(rx_bits, imageSize, bitsPerPixel);

% Plot images
figure
subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
disp(BER);



