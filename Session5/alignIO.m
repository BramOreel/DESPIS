function [out_aligned] = alignIO(out,fs)
    pulse = sin(2*pi*400*(0:1/fs:11/fs))';
    [xc,lags] = xcorr(pulse,out);
    [~,t] = max(abs(xc)); %t = de tijdsvertraging
    out_aligned = out(length(pulse)+(t+1)-50:end);
end