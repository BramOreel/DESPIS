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
T = 10;

%% Construct QAM symbol stream.
% Data blocks
[bitStream, imageData, colorMap, imageSize, bitsPerPixel] = imagetobitstream('image.bmp');
qamStream = qam_mod(bitStream, M);
streamLength = length(bitStream);



% Bitloading
ON_OFF_mask = ones(N/2-1,1); % Default all bins to one for regular transmission
bitloading_flag = 0; % If 1 then on-off/adaptive bitloading is enabled.
bitloading_type = "adaptive"; % on-off or adaptive 
BWusage = 60; % Fraction of bins to use for on-off bitloading
Nswitch = (Lt+Ld)*(N+Lcp); % The simulated channel changes every Nswitch number of samples.
smoothing_factor = .99; % Smoothing factor for simulated channel (see simulate_channel.m)

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
    H = [0;CHANNELS ;0; flipud(conj(CHANNELS))];
    if bitloading_type == "on-off"
        % Only keep bins of high energy
        
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
       H_abs = abs(H).^2;

       b_mat = floor(log2(1 + (N/2-1)*db2pow(SNR).*H_abs()/(T*sum(H_abs(1:N/2-1),"all"))))';
       bit_column = sum(b_mat(1:N/2-1),"all");
       M_vary = 2.^b_mat;

       %Rewrite to qamStream to account for these different M sizes
       j = 1; %Where are we in the bitstream?
        qamStream = [];
        while j <= streamLength
            for Np = b_mat(1:N/2-1)
            M = 2^Np;
        
            if(M == 1) && j < streamLength
                qamStream = [qamStream ; 0];
        
        
            elseif M ~=1 && j + Np - 1 <= streamLength
                 qamStream = [qamStream ; qam_mod(bitStream(j:j+Np-1),M)]; %vreemde N HIER
                j = j+Np;
            elseif M ~= 1 && j <= streamLength
                qamStream = [qamStream ; qam_mod(bitStream(j:end),M)];
                j = j+Np;
            end
            end
        end
        
        
        padLength = abs(mod(size(qamStream,1),N/2-1) -(N/2-1)) ; % Number of bits to append such that it can be divided nicely into the M-ary QAM format
        if(padLength == (N/2-1))
            padLength = 0;
        end
        qamStream = [qamStream ;zeros(padLength,1)];

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

if bitloading_type == "adaptive" && bitloading_flag == 1

    %convert the data back to a bitstream

j  = 1;
bitMatrix = [];
while j <= length(qamStream)
    for Np = b_mat(1:N/2-1)
    M = 2^Np;
    if M ~=1 && j == length(qamStream)
        nb_lastBits = streamLength - length(bitMatrix);
        lastBits_padded = qam_demod(rx_qam(j), M, Np);
        bitMatrix = [bitMatrix; lastBits_padded(end - nb_lastBits + 1:end)];
    elseif M ~= 1 && j < length(qamStream)
        bitMatrix = [bitMatrix ; qam_demod(rx_qam(j),M,Np)];
    end
    j = j+1;

    end
end
rx_bits = bitMatrix(1:streamLength);
else

    rx_bits = qam_demod(rx_qam,M,streamLength);
end





%% Bit error rate
BER = ber(rx_bits,bitStream );
disp("BER: " + BER)

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
        %xlim([0 150]);
        
        
        subplot(2,2,3);
        plot(abs(freq_res_MASK).^2); % plot the channel frequency response
        xlabel('')
        ylabel('')

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







