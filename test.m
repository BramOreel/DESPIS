% Given: channel impulse response h
h = load("Session4\channel_session4.mat").h';
N = length(h);                % Length of the impulse response
H = fft(h, N);                % Frequency response of the channel
H_abs = abs(H);               % Magnitude of the frequency response
Fs = 16000;
f = (0:N-1) * (Fs / N);       % Frequency vector, assuming sampling frequency Fs

% Plot to visualize attenuation
figure;
plot(f, 20*log10(H_abs));     % Plot magnitude in dB
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Channel Frequency Response');

%in the test function low frequencies smaller then 500 are bad
%4600 to 5760 are bad
%7466 to 8533 are bad
%10240 to 11306 are bad
%above 14750 is bad

