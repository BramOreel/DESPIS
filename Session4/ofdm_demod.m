function [ data_seq, CHANNELS ] = ofdm_demod( OFDM_seq, N, Lcp, varargin)
% OFDM demodulation
%
% INPUT:
% OFDM_seq      T1X1            Time domain OFDM sequence of length T1 samples.
% N             1X1             Total number of symbols in a single OFDM frame.
% Lcp           1X1             Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% varargin  
% For Session 3 
% empty
%
% For Session 4
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% channel       T2X1            Impulse response of channel.
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
%                               (you can ignore this until exercise 4.3)
% equalization  1X1             If 1 channel equalization is performed, if 0 no
%                               channel equalization is performed.
%
% For session 5
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
%                               (you can ignore this until exercise 4.3)
% trainblock    T3X1            Training block of T3 QAM symbols (you can
%                               ignore this until exercise 6.1.1)
%
% For Session 6
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
%                               (you can ignore this until exercise 4.3)
% trainblock    T3X1            Training block of T3 QAM symbols (you can
%                               ignore this until exercise 6.1.1)
% Lt            1X1             Number of training frames (you can ignore this
%                               until exercise 6.1.1)
% Ld            1X1             Number of data frames (you can ignore this
%                               until exercise 6.1.1)
% nbPackets     1X1             Number of packets, where one packet consist of
%                               a training and data frame. (you can ignore 
%                               this until exercise 6.1.1)
% Session 7
% Lt            1X1             Number of training frames
% M             1X1             QAM-ary constellation size
% trainblock    T3X1            Training block of T3 QAM symbols
% ON_OFF_mask   (N/2-1)X1       Mask denoting the bins to use in range (X_1...X_(N/2-1)) with a 1 denothing
%                               that a bin should be used and a 0 denoting that a bin should not be used.
% mu            1X1             Adaptive filter stepsize
% alpha         1X1             Adaptive filter regularisation factor
% type          1X1             Adaptive filter type (supported: 'NLMS')
%
% OUTPUT:
% data_seq      T4X1            QAM sequence of T4 symbols.
% CHANNELS      N/2-1XP         Frequency domain estiamted channel with P either nbPackets or 1. (
%                               you can ignore this until exercise 5.1.4)    

%% Extract input arguments
if nargin == 7 % Session 4
    streamLength = varargin{1};
    channel = varargin{2};
    ON_OFF_mask = varargin{3};
    equalization = varargin{4};
elseif nargin == 6 % Session 5
    streamLength = varargin{1};
    ON_OFF_mask = varargin{2};
    trainblock = varargin{3};
elseif nargin == 9 % Session 6
    streamLength = varargin{1};
    ON_OFF_mask = varargin{2};
    trainblock = varargin{3};
    Lt = varargin{4};
    Ld = varargin{5};
    nbPackets = varargin{6};
elseif nargin == 10 % Session 7
    Lt = varargin{1};
    M = varargin{2};
    trainblock = varargin{3};
    ON_OFF_mask = varargin{4};
    mu = varargin{5};
    alpha = varargin{6};
    type = varargin{7};
else
    equalization = 0;
end

%% Perform OFDM demodulation
% Reshape the received OFDM sequence (serial to parallel conversion)

padLength = abs(mod(size(OFDM_seq,1),N+Lcp) - (N+Lcp)) ;
    if(padLength == N+Lcp)
        padLength = 0;
    end
    OFDM_seq = [OFDM_seq; zeros(padLength,1)];
OFDM_matrix = reshape(OFDM_seq,N+Lcp,[]); %N+LcpXP

% Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
OFDM_matrix = OFDM_matrix(Lcp+1:end, :); %NxP

% Apply fft operation
QAM_matrix = fft(OFDM_matrix,N);

% Apply channel equalisation (you can ignore this until exercise 4.2.3)
CHANNELS = 1;

if equalization
    p = size(OFDM_matrix,2);
    pad = abs(mod(length(channel),N*p)-N*p);
    if pad == N*p
        pad =0;
    end
    channel = [channel; zeros(pad,1)];
    channel = reshape(channel,N,[]);
    CHANNEL = fft(channel); 
    QAM_matrix = QAM_matrix./CHANNEL;
    QAM_matrix(isinf(QAM_matrix)) = 0;
    QAM_matrix(isnan(QAM_matrix)) = 0;
end

%{
if equalization
    CHANNEL = fft(channel,N); 
    QAM_matrix = QAM_matrix./CHANNEL;
    QAM_matrix(isinf(QAM_matrix)) = 0;
    QAM_matrix(isnan(QAM_matrix)) = 0;
end
%}

%{
    padLength = abs(mod(size(channel,1),length(OFDM_seq))-(N)) ;
    if(padLength == N)
        padLength = 0;
    end
    channel = [channel; zeros(padLength,1)];
    CHANNEL = fft(channel);
    CHANNELS = reshape(CHANNEL,N,[]);
    for i = 1: size(QAM_matrix,2)
        QAM_matrix(:,1) = QAM_matrix(:,1)./CHANNELS(:,1);
    end
    QAM_matrix(isinf(QAM_matrix)) = 0;
    QAM_matrix(isnan(QAM_matrix)) = 0; %NxP
end 


if equalization
    padLength = abs(mod(size(channel,1),N*size(OFDM_matrix,2)) - N*size(OFDM_matrix,2)) ;
    if(padLength == N*size(OFDM_matrix,2))
        padLength = 0;
    end
    h = [channel; zeros(padLength,1)];
    h = reshape(h,N,[]);%NxP
    CHANNELS = fft(h);
    for i = 1: size(QAM_matrix,2) %voor elke frequentiebin P
        QAM_matrix(:,i) = QAM_matrix(:,i)./CHANNELS(:,i); %NxP
        QAM_matrix(isinf(QAM_matrix)) = 0;
        QAM_matrix(isnan(QAM_matrix)) = 0; %NxP
    end
    CHANNELS = CHANNELS(2:(N/2),:);
end
%}
    

% Remove the redundant parts of QAM_matrix 
QAM_matrix = QAM_matrix(2:(N/2),:);   %this is the fft output and needs to be scaled with an inverse %N/2-1xP


% Apply on-off mask (you can ignore this until exercise 4.3)

%QAM_matrix = QAM_matrix.*ON_OFF_mask;

% Supply streamLength number of symbols (you can ignore this until exercise 4.2)
QAM_matrix = QAM_matrix(:);
data_seq = QAM_matrix(1:streamLength);


epsilon = 10^-10;  % Smallest positive normalized number

% Loop backwards through the array to find the last valid element
idx = length(data_seq);  % Start at the last element
while idx > 0 && (abs(real(data_seq(idx))) < epsilon && abs(imag(data_seq(idx))) < epsilon)
    idx = idx - 1;  % Move backwards to the previous element
end

% Truncate the array up to the last valid element
data_seq = data_seq(1:idx);

end