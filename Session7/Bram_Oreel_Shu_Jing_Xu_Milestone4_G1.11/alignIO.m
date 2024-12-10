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


    T = 1/fs;
    f0 =  100;
    f1 = 2000;
    duration = 1; %Duration of the pulse
    t = 0:T:duration- T;
    pulse = chirp(t,f0,duration,f1)'; %1kHz sin
    sync_pulse = [pulse; zeros(fs,1)];


%% Align I/O
safety_margin = 100; % Safety margin [samples]

[correlation,lags] = xcorr(out,sync_pulse,'none'); % Find cross-correlation between output and synchronisation pulse
[~,maxIdx] = max(abs(correlation)); % Find index of maximum cross-correlation
delay = lags(maxIdx); % Find delay corresponding to that lag of maximum cross-correlation

startIdx = delay + length(sync_pulse) - safety_margin; % Find start index
out_aligned = out(startIdx:end);
%plot(out_aligned);

end
