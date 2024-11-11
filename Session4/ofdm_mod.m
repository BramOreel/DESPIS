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
    ON_OFF_mask = ones(N,1);
end

%% Construct the OFDM sequence

% Put the QAM symbols into matrix of N/2-1 rows

%display(length(QAM_seq), 'length QAM_seq before padding')
padLength = abs(mod(size(QAM_seq,1),(N/2)-1) -((N/2)-1)) ; % Number of bits to append such that it can be divided nicely into the M-ary QAM format
if(padLength == N/2-1)
    padLength = 0;
end
QAM_seq = [QAM_seq; zeros(padLength,1)];

%display(length(QAM_seq), 'length QAM_seq after padding')

QAM_matrix = reshape(QAM_seq,N/2-1,[]); %15x1280 met N = 32, N/2-1 = 15


% Construct the OFDM frames according to Figure 2 in session 3
fOFDM_frame = [zeros(1,size(QAM_matrix,2)) ; QAM_matrix ; zeros(1,size(QAM_matrix,2)) ; conj(flipud(QAM_matrix)) ];
fOFDM_frame = fOFDM_frame.*ON_OFF_mask;
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

