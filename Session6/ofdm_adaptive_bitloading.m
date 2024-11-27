function[bitMatrix] = ofdm_adaptive_bitloading(bitStream,N, Lcp, channel,SNR, b_mat,MASK, Lt, Ld,train_block)



streamLength = length(bitStream);
%MASK = ones(1,N/2-1);


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

%{
M = 2^5;
Nq = log2(M);
padLength = abs(mod(size(bitStream,1),Nq) -Nq) ;
if(padLength == Nq)
    padLength = 0;
end
bitStream = [bitStream; zeros(padLength,1)];
qamStream = qam_mod(bitStream,M);
%}


padLength = abs(mod(size(qamStream,1),N/2-1) -(N/2-1)) ;
if(padLength == (N/2-1))
    padLength = 0;
end
qamStream_padded = [qamStream ;zeros(padLength,1)];


[ofdmStream,nbPackets] = ofdm_mod(qamStream_padded,N,Lcp, MASK, Lt, Ld, train_block);


%Send the ofdm stream through the channel
rxOfdmStream = fftfilt(channel,ofdmStream);

rxOfdmStream = awgn(rxOfdmStream,SNR,'measured');
%ofdm demodulation

[rxQamStream, CHANNELS] = ofdm_demod(rxOfdmStream, N, Lcp, length(qamStream), MASK, train_block,Lt,Ld, nbPackets);



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


%bitMatrix = qam_demod(rxQamStream,M,length(bitStream));

end