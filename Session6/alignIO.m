function [ out_aligned ] = alignIO( out, fs )
% Aligns output based on a syncronization pulse.
%
% INPUT:
% out           T1x1    Recorded stream of T1 samples.
% fs            1X1     Sampling frequency [Hz].
% 
% OUTPUT:
% out_aligned   T2X1    Aligned output of T2 samples.

%% Define syncronization pulse
sync_pulse = ;

%% Align I/O
safety_margin = 50; % Safety margin [samples]

[] = xcorr(); % Find cross-correlation between output and synchronisation pulse
[~,maxIdx] = max(); % Find index of maximum cross-correlation
delay = ; % Find delay corresponding to that lag of maximum cross-correlation
startIdx = ; % Find start index
out_aligned = out(); % Align output
end

