function [BS_filtered_signal] = BS_white_noise(duration, fs, sig)


N = duration*fs; %dit is ook de filter order in fir1

%lower cutoff frequency w1 and higher cutoff frequency w2
w1 = 700/fs;
w2 = 3000/fs;

filter = fir1(N,[w1 w2],'stop');

%filterdft = abs(fft(filter));%maar de helft nodig
%filterdft = filterdft(1:floor(N/2)+1)

%plot(0:2*fs/N:fs,filterdft)

tot = fftfilt(filter,sig); %tijdsdomein


freq_totdft = fft(tot);
freq_totdft = freq_totdft(1:N/2+1);
magn = abs(freq_totdft);

figure;
title('filtered signal')
plot((0:2*fs/N:fs),mag2db(magn))
ylabel('magnitude [dB]')
xlabel('frequency [Hz]')

BS_filtered_signal = tot;
end