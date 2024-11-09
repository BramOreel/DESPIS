function [ bit_seq ] = qam_demod( QAM_seq, M, streamLength,x)
% Demodulates M-ary QAM symbols to bits. 
%
% INPUT:
% QAM_seq       T2X1    Modulated bit sequence into M-aray QAM format of length
%                       T2 symbols.
% M             1X1     M-aray QAM format (corresponding to an integer power of 2)
% varargin:
% streamLength  1X1     Length of returned bit sequence [samples].
%
% OUTPUT:
% bit_seq       T1x1    Bit sequence of T1 bits 

%% Check M is an integer power of 2
assert(sum(nextpow2(M)==log2(M))==length(M),'M is not a power of 2.')

%% Demoludation by calling qamdemod
QAM_seq = QAM_seq.*x;

bit_seq = qamdemod(QAM_seq,M);

binaryStrings = dec2bin(bit_seq'); % Convert to binary strings
% newbinaryStrings = zeros(size(binaryStrings,1),streamLength/size(binaryStrings,1));
% %add zeros so the string works

% wntdlen = streamLength/size(binaryStrings,1);
% newBin = zeros(size(binaryStrings,1),streamLength/size(binaryStrings,1));
%  for i = 1:size(binaryStrings,1)
%      actlen = size(binaryStrings(i,:),2);
%      str = str2double(cellstr(binaryStrings(i,:)))';
% 
%     if actlen < wntdlen
%         newBin(i,:) = [zeros(1,wntdlen-actlen), binaryStrings(i,:)];
% 
% 
%     end
%  end
% % 
 binaryArray = reshape(binaryStrings', [], 1); % Reshape to a single column




binaryArray = str2double(cellstr(binaryArray));




bit_seq = binaryArray(1:streamLength);

end

