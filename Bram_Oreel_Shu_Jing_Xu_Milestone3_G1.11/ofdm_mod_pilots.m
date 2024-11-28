function [OFDM_seq, nbOFDMsymb] = ofdm_mod_pilots(QAM_seq, N, Lcp, trainblock)
% OFDM modulation baed on pilots.
%
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% trainblock    T2X1        Training block of T2 QAM symbols 
%
% OUTPUT:
% data_seq      T3X1        QAM sequence of T3 symbols.
% nbOFDMsymb    1X1         Number of OFDM frames.


carriers_used = N/4-1;



%display(length(QAM_seq), 'length QAM_seq before padding')
padLength = abs(mod(size(QAM_seq,1),carriers_used) - carriers_used) ; % Number of bits to append such that it can be divided nicely into the M-ary QAM format
if(padLength == carriers_used)
    padLength = 0;
end
QAM_seq = [QAM_seq; zeros(padLength,1)];

%Place training data on the odd bins i guess, really dont wanna do this,
%make a full vector and just zero out the frequencies you dont want
mat = [1 ; 0];
mask = repmat(mat,ceil((N/2-1)/2),1);
mask = mask(1:end-1);

QAM_matrix = [];
i = 1;
%Place data in the zero tones until every data piece has been filled
while i <= length(QAM_seq)
       tonevec = trainblock.*mask;
       for j = 1:N/2-1
            if tonevec(j) == 0
                tonevec(j) = QAM_seq(i);
                i = i +1;
            end
       end
       QAM_matrix = [QAM_matrix , tonevec];
end






%% Extract data
%QAM_matrix = ; % even bins for data    

nbOFDMsymb = size(QAM_matrix,2);
% Construct the OFDM frames
fOFDM_frame = [zeros(1,size(QAM_matrix,2)) ; QAM_matrix ; zeros(1,size(QAM_matrix,2)) ; conj(flipud(QAM_matrix)) ];

% Apply the inverse Fourier transform (IFFT)
OFDM_frame =  ifft(fOFDM_frame);

% Add in the cyclic prefix
OFDM_frame = [ OFDM_frame(end-Lcp+1:end, :) ;OFDM_frame] ;

% Serialize the set of OFDM frames
OFDM_seq = OFDM_frame(:);
end

