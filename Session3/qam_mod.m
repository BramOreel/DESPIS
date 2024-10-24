function [ QAM_seq ] = qam_mod( bit_seq, M )
% Modulates a bit sequence into M-ary QAM format. 
%
% INPUT:
% bit_seq   T1x1    Bit sequence of T1 bits 
% M         1X1     M-aray QAM format (corresponding to an integer power of 2)
%
% OUTPUT:
% QAM_seq   T2X1    Modulated bit sequence into M-aray QAM format of length
%                   T2 symbols.

%% Check M is an integer power of 2
%assert(sum(nextpow2(M)==log2(M))==length(M),'M is not a power of 2.')

%% Append zero-bits to the end of the sequence to make the sequence 
N = log2(M); % Number of bits per QAM symbol
padLength = abs(mod(size(bit_seq,1),N) -N) ; % Number of bits to append such that it can be divided nicely into the M-ary QAM format
bit_seq = [bit_seq; zeros(padLength,1)]; % Padded bit sequence

% Check that the bit sequence is correctly padded
assert( mod(length(bit_seq),sum(N))==0,'Bit sequence should contain a number of bits that is a multiple of log2(M).')

%% Call to qammod() to obtain the M-ary QAM symbols
QAM_seq = qammod(bit2int(bit_seq,N),M);
x = sum(abs(QAM_seq),"all")/size(QAM_seq,1)
QAM_seq = QAM_seq./x;

end

