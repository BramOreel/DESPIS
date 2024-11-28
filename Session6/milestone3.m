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

max_Nq = 4; % Maximum QAM size (64-QAM)

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);
streamLength = length(bitStream);



% Bitloading
ON_OFF_mask = ones(1,N/2-1); % Default all bins to one for regular transmission
bitloading_flag = 1; % If 1 then on-off/adaptive bitloading is enabled.
bitloading_type = "on-off"; % on-off or adaptive 
BWusage = 60; % Fraction of bins to use for on-off bitloading
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

%% Determine bit loading
% Channel estimation based on a dummy transmission for bitloading

train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits
train_block = qam_mod(train_bits,M); % QAM modulate -> (N/2-1) rijen

%if the bitloading flag in toggled, we will perform two transmissions
%Basically just the train block about 10 times (fs*2/N blocks)
%[trainStream, nbPackets] = ofdm_mod(train_block,N,Lcp,ON_OFF_mask, floor(fs*2/N),0,train_block); %probleem: Ld = 0 geeft nbPackets = Inf
[trainStream, nbPackets] = ofdm_mod(train_block,N,Lcp,ON_OFF_mask, Lt,Ld,train_block);

if bitloading_flag

    % Dummy transmission
    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session6.mat',smoothing_factor);
        aligned_Rx = awgn(aligned_Rx,SNR,'measured'); %y = h*x+n -> bestaat uit trainingsdata en dummy data
    else
        
        duration = 1;
        T = 1/fs;
        t = 0:T:duration- T;
        pulse = chirp(t,100,duration,2000)';
        sync_pulse = [pulse; zeros(fs,1)];

        [simin,nbsecs] = initparams(trainStream,fs,sync_pulse);
        size(simin)
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO(Rx,fs);
    end
    [aligned_RX, CHANNELS] = ofdm_demod(aligned_Rx,N,Lcp, length(train_block),ones(1,N/2-1),train_block,Lt,Ld,nbPackets);

    if bitloading_type == "on-off"
        % Only keep bins of high energy
        H = [0;CHANNELS ;0; flipud(conj(CHANNELS))];
        H_abs = abs(H);

        [H_sorted, idx_sorted] = sort(H_abs(1:N), 'descend');

                                % Percentage of frequency bins to use
        num_bins = floor(length(idx_sorted) * (BWusage / 100));        % Number of bins to use
        selected_bins = idx_sorted(1:num_bins);      % Indices of selected bins, these bins have a good SNR

        % Create a frequency mask to use only the selected bins
        frequency_mask = zeros(N, 1);
        frequency_mask(selected_bins) = 1;
        frequency_mask = frequency_mask(1:N/2-1);
        ON_OFF_mask = frequency_mask; % ON-OFF mask with 1 denoting the usage of a bin.

    elseif bitloading_type == "adaptive"
        H = CHANNELS;
        %SNR in het frequentiedomein
        
        verschil = abs(length(H)-length(train_block));
        if length(train_block) < length(H)
            train_block = [train_block; zeros(verschil,size(train_block,2))];
            aligned_RX = [aligned_RX; zeros(verschil,size(aligned_RX,2))]; %padden met zeros als dummy_qam niet lang genoeg is, maar dit geeft meer BER!
        elseif length(train_block) > length(H)
            H = [H;zeros(verschil,size(H,2))];
        end
        
        NOISE = aligned_RX - H.*train_block; 
        PSDn = (abs(NOISE).^2)/(N*fs);

        T = 10;
        
        SNR_per_bin = (abs(H).^2)./(T*PSDn);

        [sorted_SNR, sorted_indices] = sort(SNR_per_bin, 'descend');
        num_active_tones = ceil(BWusage / 100 * (N/2-1));
        used_indices = sorted_indices(1:num_active_tones);

        active_tones = zeros(N/2-1, 1);
        active_tones(used_indices) = 1;
        
        %Shannon
        b_mat = floor(log2(1+ SNR_per_bin)); 
        for i=1:length(active_tones)
            if b_mat(i)>max_Nq
                b_mat(i) = max_Nq;
            end
        end

        M_vary = 2.^b_mat;     % Constellation sizes
        
        ch = load('channel_session6.mat').h';
        
        Rx_bitstream = ofdm_adaptive_bitloading(bitStream,N, Lcp,ch,SNR, b_mat.',active_tones, Lt, Ld, train_block);

        
        BER = ber(Rx_bitstream,bitStream );

        % Construct image from bitstream
        imageRx = bitstreamtoimage(Rx_bitstream, imageSize, bitsPerPixel);

        % Plot images
        figure
        subplot(2,1,1); colormap(colorMap); image(imageData); axis image; title('Original image'); drawnow;
        subplot(2,1,2); colormap(colorMap); image(imageRx); axis image; title('Received image'); drawnow;
        disp(BER);

        %Noise bepalen
        %M bepalen
        
        

    end
end




%% OFDM modulation
if bitloading_type == "on-off"
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
scatterplot(rx_qam);
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
end




%% Visualize_demod

% Transmission of a picture across the channel


%% Construct QAM symbol stream.


%% Determine bit loading
% Channel estimation based on a dummy transmission for bitloading

train_bits = randi([0 1],Nq*(N/2-1),1); % Generate a random vector of bits
train_block = qam_mod(train_bits,M); % QAM modulate

%if the bitloading flag in toggled, we will perform two transmissions
%Basically just the train block about 10 times (fs*2/N blocks)
trainStream = ofdm_mod(train_block,N,Lcp,ON_OFF_mask, floor(fs*2/N),0,train_block);


if bitloading_flag

    % Dummy transmission
    if channel == "simulation"
        aligned_Rx = simulate_channel(trainStream, Nswitch,'channel_session6.mat',smoothing_factor);
        aligned_Rx = awgn(aligned_Rx,SNR,'measured');
    else
        
        duration = 1;
        T = 1/fs;
        t = 0:T:duration- T;
        pulse = chirp(t,100,duration,2000)';
        sync_pulse = [pulse; zeros(fs,1)];

        [simin,nbsecs] = initparams(trainStream,fs,sync_pulse);
        size(simin)
        sim('recplay');
        Rx = simout.signals.values(:,1);
        aligned_Rx = alignIO(Rx,fs);
    end

    % Extract the estimated channels
    [~, CHANNELS] = ofdm_demod(aligned_Rx,N,Lcp, 0,ones(1,N/2-1),train_block,floor(fs*2/N),0,1);

    if bitloading_type == "on-off"
        % Only keep bins of high energy
        H = [0;CHANNELS ;0; flipud(conj(CHANNELS))];
        H_abs = abs(H);
        [H_sorted, idx_sorted] = sort(H_abs(1:N), 'descend');

                                % Percentage of frequency bins to use
        num_bins = floor(length(idx_sorted) * (BWusage / 100));        % Number of bins to use
        selected_bins = idx_sorted(1:num_bins);      % Indices of selected bins, these bins have a good SNR

        % Create a frequency mask to use only the selected bins
        frequency_mask = zeros(N, 1);
        frequency_mask(selected_bins) = 1;
        frequency_mask = frequency_mask(1:N/2-1);
        ON_OFF_mask = frequency_mask; % ON-OFF mask with 1 denoting the usage of a bin.
    elseif bitloading_type == "adaptive"
       % M = ;     % Constellation sizes
    end
end



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

%Data rate = N*Nq*R with R = 1/Tlcp + Tdata
%It takes Rx / fs  seconds for all packets
%This divided by nbpackets is our refresh rate

t_packet = length(aligned_Rx)/(fs*nbPackets); %for LD + Lt packets; % Duration of one packet [s]

idx = 1;
while idx <= nbPackets
    try
        CHANNEL = CHANNELS(:,idx);
        CHANNEL_MASK = CHANNEL.*ON_OFF_mask;
        freq_res = [0;CHANNEL ;0; flipud(conj(CHANNEL))];
        freq_res_MASK = [0;CHANNEL_MASK ;0; flipud(conj(CHANNEL_MASK))];
        h_est = ifft(freq_res,N);

        

       
     

        figure(1);
        
        subplot(2,2,1);
        plot(h_est);
        ylim([-0.05 0.05]); % plot the channel impulse response
        xlim([10 150]);
        xlabel('')
        ylabel('')
        title("channel impulse respons")
        xlabel('Time [samples]')
        ylabel('Impulse response [arb.]')
        %xlim([0 150]);
        
        
        subplot(2,2,3);
        plot(abs(freq_res_MASK).^2); % plot the channel frequency response
        xlabel('')
        ylabel('')
        ylim([0 0.03])
        title("channel frequency response")
        xlabel('Frequency [Hz]')
        ylabel('Magnitude response [dB]')
        subplot(2,2,2);
        colormap(colorMap); image(imageData); axis image; title('Original image'); % the original image

        
        subplot(2,2,4);
        demoded_bitstream = rx_bits(1:Nq*(N/2-1)*Ld*idx); % the bit stream that is received by the receiver so far
        imageRx = bitstreamtoimage(demoded_bitstream, imageSize, bitsPerPixel);
        colormap(colorMap); image(imageRx); axis image; title(['Received image']); drawnow;
        drawnow;
    catch
        break;
    end
    
    pause(t_packet);
    idx = idx + 1;
end

% 
% A static or slowly changing impulse response and frequency response indicate a stable channel, suitable for high-speed communication.
% Rapid time variations suggest a challenging environment, such as a mobile user or dynamic obstacles, requiring robust channel estimation and equalization techniques.
% A frequency response with deep fading notches implies that certain subcarriers may experience severe attenuation, necessitating error correction or adaptive modulation schemes.



