
clear all;close all;
L = 1025;
SNR = 10;
bit_seq = randi([0, 1], 1,L)' ;
v = qam_mod(bit_seq,64);
req_QAM_seq  = awgn(v,SNR);  

w = qam_demod(v,64,L);




%scatterplot(req_QAM_seq)