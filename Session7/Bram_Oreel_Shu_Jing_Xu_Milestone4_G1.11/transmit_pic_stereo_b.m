% Transmission of a picture across the channel
%% Cleanup
clear; close all; clc;

%% Parameters.
N = 1024; % Total number of symbols in a single OFDM frame, i.e., the DFT size
Lcp = 300; % Cyclic prefix length [samples].
M = 8; % QAM constellation size.
Nq = log2(M);
SNR = 15; % SNR of transmission [dB].
Lt = 10; % Number of training frames.
Ld = 10; % Number of data frames.
fs = 16000; % Sampling frequency [Hz].
channel = "acoustic"; % simulation or acoustic
Equalization = "packet"; % Equalization mode (see ofdm_demod_stereo.m)
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)
mu = 0.02;
alpha = 1;

%% Initial channel estimation.
% Construct train block.
train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits
train_frame = qam_mod(train_bits,M); % QAM training symbols

% Initial channel estimation.
ofdm_train_seq = ofdm_mod([],N,Lcp,Lt,train_frame,ones(N/2-1,1));

% Construct training sequence for transmitter side channel estimation.

ofdm_train_seq_stereo = [ofdm_train_seq zeros(length(ofdm_train_seq),1); zeros(length(ofdm_train_seq),1) ofdm_train_seq];

%% Transmit OFDM sequence.
if channel == "simulation"
    [simin,nbsecs] = initparams_stereo(ofdm_train_seq_stereo,fs);
    Rx = simulate_channel_stereo(simin, Nswitch,'channel_stereo_session7.mat',smoothing_factor);
    Rx = awgn(Rx,SNR,'measured');
    aligned_Rx = alignIO(Rx,fs);
elseif channel == "acoustic"
    [simin,nbsecs] = initparams_stereo(ofdm_train_seq_stereo,fs);
    sim('recplay');
    Rx = simout.signals.values(:,1);
    aligned_Rx = alignIO(Rx,fs);
end

timesig_left = aligned_Rx(1:length(ofdm_train_seq));
timesig_right = aligned_Rx(length(ofdm_train_seq)+1: 2*length(ofdm_train_seq));

[~,channel_left_freq] = ofdm_demod_stereo(timesig_left,N,Lcp,train_frame,Lt,0,M,0,"fixed",mu,alpha);
[~,channel_right_freq] = ofdm_demod_stereo(timesig_right,N,Lcp,train_frame,Lt,0,M,0,"fixed",mu,alpha);

freq_res_left = [0;channel_left_freq ;0; flipud(conj(channel_left_freq))];
channel_left_time = ifft(freq_res_left,N);

freq_res_right = [0;channel_right_freq ;0; flipud(conj(channel_right_freq))];
channel_right_time = ifft(freq_res_right,N);

H = [channel_left_freq(1:N/2-1),channel_right_freq(1:N/2-1)];


% Plot channel estimations.
figure(1);
subplot(2,1,1);
plot(channel_left_time(1:Lcp));
hold on
plot(channel_right_time(1:Lcp));
xlabel('Samples')
ylabel('Magnitude [arb.]')   
title('Impulse response.');
legend('left','right');
subplot(2,1,2);
plot(abs(freq_res_left).^2);
hold on
plot(abs(freq_res_right).^2);
title('Frequency response.');
xlabel('Samples')
ylabel('Magnitude [dB]')   
legend('left','right');

%% Transmission, including beamforming.
% Construct QAM symbol stream.
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M); %Our qam sequence we will use during the adaptive filter
streamLength = length(bitStream);

% OFDM modulation.
[Tx, a, b, nbPackets] = ofdm_mod_stereo(qamStream,N,Lcp,H,Lt,Ld,train_frame,"packet");

channel = "simulation";

% Transmit data sequence.
if channel == "simulation"
    [simin,nbsecs] = initparams_stereo(Tx,fs);
    Rx = simulate_channel_stereo(simin, Nswitch,'channel_stereo_session7.mat',smoothing_factor);
    Rx = awgn(Rx,SNR,'measured');
    aligned_Rx = alignIO(Rx,fs);
    aligned_Rx = aligned_Rx(1:length(Tx)).';
end

%% Post transmission processing
% OFDM demodulation.
[qam_seq, CHANNELS] = ofdm_demod_stereo(aligned_Rx,N,Lcp,train_frame,Lt,Ld,M,nbPackets,"packet",mu,alpha);
scatterplot(qam_seq)
CHANNEL_combo = CHANNELS(:,1); % Extract first estimated channel
figure
plot(abs(CHANNEL_combo).^2)

% QAM demodulation.
rx_bits = qam_demod(qam_seq,M,streamLength);

% Display results Construct image from bitstream
figure;
imageRx = bitstreamtoimage(rx_bits, imageSize, bitsPerPixel);
% Plot images
subplot(1,2,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
subplot(1,2,2); colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
% Calculate BER
BER = ber(rx_bits,bitStream) ;
%}