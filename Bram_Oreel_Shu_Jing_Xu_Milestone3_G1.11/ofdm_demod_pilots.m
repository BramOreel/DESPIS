function [data_seq, CHANNELS] = ofdm_demod_pilots( OFDM_seq, N, Lcp, streamLength, trainblock,nbOFDMsymb)
% OFDM demodulation using pilot tones
% INPUT:
% OFDM_seq      T1X1            Time domain OFDM sequence of length T1 samples.
% N             1X1             Total number of symbols in a single OFDM frame.
% Lcp           1X1             Cyclic prefix length [samples] (you can ignore this until exercise 3.2.4)
% streamLength  1X1             Length of QAM sequence after call to qammod [samples].
% varargin  
% 
% Session 5
% trainblock    T2X1            Training block of T2 QAM symbols (you can
%                               ignore this until exercise 6.1.1)
% nbOFDMsymb    1X1             Number of OFDM frames.
%
% OUTPUT:
% data_seq      T3X1            QAM sequence of T3 symbols.
% CHANNELS      N/2-1XnbOFDMsymbFrequency domain estimated channel for each frame. (
%                               you can ignore this until exercise 5.1.4)   

%% Perform OFDM demodulation
% Reshape the received OFDM sequence (serial to parallel conversion)
OFDM_matrix = reshape(OFDM_seq,N+Lcp,[]);
% Remove the cyclic prefix (you can ignore this until exercise 3.2.4)
OFDM_matrix = OFDM_matrix(Lcp+1:end, :);
% Apply fft operation
QAM_matrix = fft(OFDM_matrix,N);
% Remove the redundant parts of QAM_matrix
QAM_matrix = QAM_matrix(2:(N/2),:);

data_matrix = zeros(N/4-1,nbOFDMsymb); % Placeholder for data
CHANNELS = zeros(N/2-1,nbOFDMsymb); % Plaecholder for channels
pilot_indices = 1:2:N/2-1;
data_indices = 2:2:N/2-1;

for pIdx = 1:nbOFDMsymb % Loop across frames
    % Extract packet.
    dataFrames = zeros(N/2-1,1);
    dataFrames(data_indices) = QAM_matrix(data_indices,pIdx); %matrix for the qam data
    pilotFrames = zeros(N/4,1);  %zero matrix for the pilot tones
    pilotFrames = QAM_matrix(pilot_indices,pIdx);
        
    % Channel estimation
    CHANNELS(:,pIdx) = ofdm_channelest_pilots(pilotFrames,trainblock,N,Lcp); % Save channel call ofdm channelest pilots here

    % Equalization of data frames
    equalised = dataFrames(:)./CHANNELS(:,pIdx);
    data_matrix(:,pIdx) = equalised(data_indices)./2;
end

data_seq = data_matrix(:);
data_seq = data_seq(1:streamLength);
end



