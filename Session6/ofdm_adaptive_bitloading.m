function[bitMatrix] = ofdm_adaptive_bitloading(trainStream, N, Lcp, M_vary)

streamLength = length(trainStream);

j = 1; %Where are we in the bitstream?
qamStream = [];
while j <= streamLength
    for M = M_vary(1:end)
        if (M == 1) && j < streamLength
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
    for M = M_vary(1:end)
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