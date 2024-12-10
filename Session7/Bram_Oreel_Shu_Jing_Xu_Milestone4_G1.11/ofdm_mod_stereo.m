function [OFDM_seq,a,b,nbPackets] = ofdm_mod_stereo( QAM_seq, N, Lcp, H, Lt, Ld, trainblock, Equalization)
% Stereo OFDM modulation
% 
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples].
% H             N/2-1X2     Two frequency domain channels for a DFT size of N.
% Lt            1X1         Number of training frames.
% Ld            1X1         Number of data frames.
% trainblock    T2X1        Training block of T2 QAM symbols
% Equalization  String      Equalization mode: If "fixed" (for transmit_pic_stereo_a), training frames are
%                           transmitted first to estimate the channel and
%                           no additional updating of the channel is
%                           performed thereafter as only data frames
%                           follow. If "packet" packet-based equalization should be
%                           used. (for transmit_pic_stereo_b)
% 
% OUTPUT:
% OFDM_seq      T3X2        Two channel time domain OFDM sequence of length T3 samples.
% a             N/2-1X1     Per-bin factor associated with the first stream.
% b             N/2-1X1     Per-bin factor associated with the second stream.
% nbPackets     1X1         Number of packets, where one packet consist of
%                           a training and data frame.

%perform an analog process, except leave the Lt and Ld shit blank for now
%pad the signal
carriers_used = N/2-1;
mod_pad = mod(size(QAM_seq,1),carriers_used);
div_pad = floor(size(QAM_seq,1)/carriers_used);
if mod_pad ~=0
    div_pad = div_pad + 1; %numbers of rows in the matrix in think
end
padLength = abs(mod(size(QAM_seq,1),carriers_used) - carriers_used) ; % Number of bits to append such that it can be divided nicely into the M-ary QAM format
if(padLength == carriers_used)
    padLength = 0;
end
QAM_seq = [QAM_seq; zeros(padLength,1)];

%No on off bit loading here
QAM_matrixON = reshape(QAM_seq,carriers_used,[]);


%Code for putting Data packets and training packets in certain orders, not
%required yet
if(Equalization == "fixed")
nbPackets = size(QAM_matrixON,2);
QAM_matrix = zeros(N/2-1,Lt);
for j = 1:Lt
    QAM_matrix(:,j) =  trainblock;  
end
QAM_matrix = [QAM_matrix, QAM_matrixON];

nbPackets = size(QAM_matrixON,2);

else if(Equalization == "packet")
    nbPackets = ceil(size(QAM_matrixON,2)/Ld);
    QAM_matrix = [];
    i = 1;
    while i <= size(QAM_matrixON,2)
        for j = 1:Lt
            QAM_matrix = [QAM_matrix , trainblock];   
        end
    
        if(i + Ld <= size(QAM_matrixON,2))
            QAM_matrix = [QAM_matrix, QAM_matrixON(:,i:i+Ld-1)];
        else
            QAM_matrix = [QAM_matrix, QAM_matrixON(:,i:end),zeros(N/2-1,nbPackets*Ld - i)];
        end
    
    
        %Account for the case if the on off mask needs to be estimated
        if(Ld == 0)
            break;
        end
        i = i + Ld;
    end

else
    QAM_matrix = QAM_matrixON;
end
end

%7.3: Calculate the a & b using the channel estimation
[a,b] = fixed_transmitter_side_beamformer(H);

%Create two QAM_matrices, one scaled with a, another scaled with b
QAM_matrixA = QAM_matrix.*a;
QAM_matrixB = QAM_matrix.*b;

%Convert to ofdm streams
fOFDMA = [zeros(1,size(QAM_matrixA,2)) ; QAM_matrixA ; zeros(1,size(QAM_matrixA,2)) ; conj(flipud(QAM_matrixA)) ];
fOFDMB = [zeros(1,size(QAM_matrixB,2)) ; QAM_matrixB ; zeros(1,size(QAM_matrixB,2)) ; conj(flipud(QAM_matrixB)) ];
%time domain transform
OFDMstreamA = ifft(fOFDMA);
OFDMstreamB = ifft(fOFDMB);
% Add in the cyclic prefix
ofdmStream1 = [ OFDMstreamA(end-Lcp+1:end, :) ;OFDMstreamA];
ofdmStream2 = [ OFDMstreamB(end-Lcp+1:end, :) ;OFDMstreamB];
% Serialize the set of OFDM frames

OFDM_seq = [ofdmStream1(:),ofdmStream2(:)];

end