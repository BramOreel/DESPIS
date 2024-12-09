function [ data_seq, CHANNELS] = ofdm_demod_stereo( OFDM_seq, N, Lcp, trainblock, Lt, Ld, M, nbPackets, Equalization, mu, alpha )
% Performs stereo OFDM demodulation
% 
% INPUT:
% OFDM_seq      T1X1            Time domain OFDM sequence of length T samples.
% N             1X1             Total number of symbols in a single OFDM frame.
% Lcp           1X1             Cyclic prefix length [samples].
% trainblock    T2X1            Training block of T2 QAM symbols 
% Lt            1X1             Number of training frames 
% Ld            1X1             Number of data frames 
% M             1X1             QAM-ary constellation size
% Lt            1X1             Number of training frames
% Ld            1X1             Number of data frames
% nbPackets     1X1             Number of packets, where one packet consist of
%                               a training and data frame.
% Equalization  String          Equalization mode: If "fixed" (for transmit_pic_stereo_a), training frames are
%                               transmitted first to estimate the channel and
%                               no additional updating of the channel is
%                               performed thereafter as only data frames
%                               follow. If "packet" packet-based equalization should be
%                               used. (for transmit_pic_stereo_b)
% mu            1X1             NLMS stepsize
% alpha         1X1             NLMS regularization factor
% 
% OUTPUT:
% data_seq      T3X1            QAM sequence of T3 symbols.
% CHANNELS      N/2-1XP         Frequency domain estimated combined channel for each 1..P. 

pad = abs(mod(length(OFDM_seq),N+Lcp)-(N+Lcp));
if pad == (N+Lcp)
    pad =0;
end
OFDM_seq = [OFDM_seq;zeros(pad,1)];
% Reshape the received OFDM sequence (serial to parallel conversion)
OFDM_matrix = reshape(OFDM_seq,N+Lcp,[]);

% Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
OFDM_matrix = OFDM_matrix(Lcp+1:end, :);

% Apply fft operation
QAM_matrix = fft(OFDM_matrix,N);

%These values we actually need
QAM_matrix_usefull = QAM_matrix(2:(N/2),:); 

if(Equalization == "fixed")
    CHANNELS = [];
    QAM_matrix = [];

    QAM_matrix_train = QAM_matrix_usefull(:,1:Lt);
    QAM_matrix_train = sum(QAM_matrix_train,2)./Lt;
    QAM_matrix_data  = QAM_matrix_usefull(:,Lt+1:end);

    CHANNEL = QAM_matrix_train(:,1)./trainblock(:,1);
    CHANNELS = [CHANNELS,CHANNEL];
    Equaliser = repmat(CHANNEL,1,size(QAM_matrix_data,2));
    QAM_matrix = [QAM_matrix, QAM_matrix_data./Equaliser];




elseif(Equalization == "packet")
        i = 1;
    CHANNELS = [];
    QAM_matrix = [];
    
    ie = nbPackets*(Lt+Ld);
    while i <= ie
        e = i+Lt-1;
        QAM_matrix_train2 = QAM_matrix_usefull(:,i:e);
        QAM_matrix_train = sum(QAM_matrix_train2,2)./Lt;
    
        QAM_matrix_data  = QAM_matrix_usefull(:,i+Lt:i+Lt+Ld-1);
        CHANNEL = QAM_matrix_train(:,1)./trainblock(:,1);
    
        Equaliser = repmat(CHANNEL,1,Ld);
    
    
        QAM_matrix = [QAM_matrix, QAM_matrix_data./Equaliser];
        CHANNELS = [CHANNELS,CHANNEL];
    
        i = i + Ld + Lt;
    end
end

QAM_matrix(isinf(QAM_matrix)) = 0;
QAM_matrix(isnan(QAM_matrix)) = 0;


% Supply streamLength number of symbols (you can ignore this until exercise 4.2)
%We need to truncate the array so that the dimensions fit again


data_seq = QAM_matrix(:);

end

