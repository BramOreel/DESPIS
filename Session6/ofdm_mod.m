function [ OFDM_seq, nbPackets ] = ofdm_mod( QAM_seq, N, Lcp, varargin)
% OFDM modulation
%
% INPUT:
% QAM_seq       T1X1        Sequence containing QAM symbols
% N             1X1         Total number of symbols in a single OFDM frame.
% Lcp           1X1         Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% varargin
% 
% Session 4 and 5
% ON_OFF_mask   (N/2-1)X1   Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                           that a bin should be used and a 0 denoting that a bin should not be used.
%                           (you can ignore this until exercise 4.3)
%
% Session 6
% ON_OFF_mask   (N/2-1)X1   Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                           that a bin should be used and a 0 denoting that a bin should not be used.
%                           (you can ignore this until exercise 4.3)
% Lt            1X1         Number of training frames (you can ignore this
%                           until exercise 6.1.1)
% Ld            1X1         Number of data frames (you can ignore this
%                           until exercise 6.1.1)
% trainblock    T2X1        Training block of T2 QAM symbols (you can
%                           ignore this until exercise 6.1.1)
% Session 7
% Lt            1X1         Number of training frames 
% trainblock    T2X1        Training block of T2 QAM symbols
% ON_OFF_mask   (N/2-1)X1   Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                           that a bin should be used and a 0 denoting that a bin should not be used.
%
% OUTPUT:
% OFDM_seq      T3X1        Time domain OFDM sequence of length T3 samples.
% nbPackets     1X1         Number of packets, where one packet consist of
%                           a training and data frame. (you can ignore 
%                           this until exercise 6.1.1)

%% Extract data
if nargin == 4
    ON_OFF_mask = varargin{1};
elseif nargin == 7
    ON_OFF_mask = varargin{1};
    Lt = varargin{2};
    Ld = varargin{3};
    trainblock = varargin{4};
elseif nargin == 6
    Lt = varargin{1};
    trainblock = varargin{2};
    ON_OFF_mask = varargin{3};
else
    ON_OFF_mask = ones(1,N/2-1);
end

%% Construct the OFDM sequence


% Put the QAM symbols into matrix of N/2-1 rows


%Map the QAM symbols onto frequencies which are allowed by the mask. This
%means that our padlength is not effective anymore i think. We can however
%extract the number of allowed carrier frequencies and calculate a new
%padding

%Final dimensions are N/2-1 * div_pad
%input is length x
%the number of columns will not change, only number of rows
carriers_used = sum(ON_OFF_mask,"all");
mod_pad = mod(size(QAM_seq,1),carriers_used);
div_pad = floor(size(QAM_seq,1)/carriers_used);
if mod_pad ~=0
    div_pad = div_pad + 1; %numbers of rows in the matrix in think
end

%display(length(QAM_seq), 'length QAM_seq before padding')
padLength = abs(mod(size(QAM_seq,1),carriers_used) - carriers_used) ; % Number of bits to append such that it can be divided nicely into the M-ary QAM format
if(padLength == carriers_used)
    padLength = 0;
end
QAM_seq = [QAM_seq; zeros(padLength,1)];

%Construct the QAM matrix with the knowledge that zeros have to be inserted
QAM_matrixON = reshape(QAM_seq,carriers_used,[]);
for k = 1:length(ON_OFF_mask)
    if ON_OFF_mask(k) == 0
        if(k == 1)
            QAM_matrixON = [zeros(1,size(QAM_matrixON,2)); QAM_matrixON];
        
        elseif(k==2)
            QAM_matrixON = [QAM_matrixON(1,:); zeros(1,size(QAM_matrixON,2)); QAM_matrixON(2:end,:)];
        
        elseif(k == length(ON_OFF_mask) -1)
            QAM_matrixON = [QAM_matrixON(1:k-2,:); zeros(1,size(QAM_matrixON,2)); QAM_matrixON(k-1,:)]; 
        
        elseif(k == length(ON_OFF_mask))
            QAM_matrixON = [QAM_matrixON;zeros(1,size(QAM_matrixON,2))];
        
        else
            QAM_matrixON = [QAM_matrixON(1:k-1,:); zeros(1,size(QAM_matrixON,2)); QAM_matrixON(k:end,:)];  %voor de tweede index wil je de eerste rij saven en de rest erna kopieren

        end
    end
end

if nargin == 7

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




% Construct the OFDM frames according to Figure 2 in session 3
fOFDM_frame = [zeros(1,size(QAM_matrix,2)) ; QAM_matrix ; zeros(1,size(QAM_matrix,2)) ; conj(flipud(QAM_matrix)) ];
%eerste rij nullen voor DC-componenten
%symmetrisch om reële tijdsignaal te bekomen


% Apply the inverse Fourier transform (IFFT)
OFDM_frame = ifft(fOFDM_frame); %size = 32x1280

% Add in the cyclic prefix
OFDM_frame = [ OFDM_frame(end-Lcp+1:end, :) ;OFDM_frame]; %size = (32+Lcp)x1280 = 48x1280
%display(size(OFDM_frame),'OFDM_frame with Lcp')

% Serialize the set of OFDM frames

OFDM_seq = OFDM_frame(:);
%display(length(OFDM_seq),'length(OFDM_seq)') % 48x1280 = 61440


elseif nargin == 4

QAM_matrix = QAM_matrixON;


% Construct the OFDM frames according to Figure 2 in session 3
fOFDM_frame = [zeros(1,size(QAM_matrix,2)) ; QAM_matrix ; zeros(1,size(QAM_matrix,2)) ; conj(flipud(QAM_matrix)) ];
%eerste rij nullen voor DC-componenten
%symmetrisch om reële tijdsignaal te bekomen


% Apply the inverse Fourier transform (IFFT)
OFDM_frame = ifft(fOFDM_frame); %size = 32x1280

% Add in the cyclic prefix
OFDM_frame = [ OFDM_frame(end-Lcp+1:end, :) ;OFDM_frame]; %size = (32+Lcp)x1280 = 48x1280
%display(size(OFDM_frame),'OFDM_frame with Lcp')

% Serialize the set of OFDM frames

OFDM_seq = OFDM_frame(:);
%display(length(OFDM_seq),'length(OFDM_seq)') % 48x1280 = 61440

    
end

end

