function[bitMatrix] = ofdm_adaptive_bitloading(bitStream,N, Lcp, channel,SNR, M_vary,MASK, Lt, Ld,train_block)

streamLength = length(bitStream);


%Moduleert bitstream met varierende Nq.
j = 1; %Where are we in the bitstream?
qamStream = [];
while j <= streamLength
    for i = 1:size(M_vary,1)
        M=M_vary(i);
        Nq = log2(M);
        if (M_vary(i) == 1) && (j < streamLength)
            qamStream = [qamStream ; 0];
        elseif M_vary(i) ~=1 && j + Nq - 1 <= streamLength
            qamStream = [qamStream ; qam_mod(bitStream(j:j+Nq-1),M_vary(i))]; %vreemde N HIER
            j = j+Nq;
        elseif M_vary(i) ~= 1 && j <= streamLength
            qamStream = [qamStream ; qam_mod(bitStream(j:end),M_vary(i))];
            j = j+Nq;
        end
    end
end



padLength = abs(mod(size(qamStream,1),N/2-1) -(N/2-1)) ;
if(padLength == (N/2-1))
    padLength = 0;
end
qamStream_padded = [qamStream ;zeros(padLength,1)];


[ofdmStream,nbPackets] = ofdm_mod(qamStream_padded,N,Lcp,MASK, Lt, Ld, train_block);


%Send the ofdm stream through the channel
rxOfdmStream = fftfilt(channel,ofdmStream);

rxOfdmStream = awgn(rxOfdmStream,SNR,'measured');
%ofdm demodulation

[rxQamStream, CHANNELS] = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), MASK, train_block,Lt,Ld, nbPackets);


%convert the data back to a bitstream
j  = 1;
bitMatrix = [];
while j <= length(rxQamStream)
    for i = 1:size(M_vary,1)
        M = M_vary(i);
        Nq = log2(M);
        if M ~=1 && j == length(rxQamStream)
            nb_lastBits = streamLength - length(bitMatrix);
            lastBits_padded = qam_demod(rxQamStream(j), M, Nq);
            bitMatrix = [bitMatrix; lastBits_padded(end - nb_lastBits + 1:end)];
        elseif M ~= 1 && j < length(rxQamStream)
            bitMatrix = [bitMatrix ; qam_demod(rxQamStream(j),M,Nq)];
        end
    j = j+1;

    end
end

end