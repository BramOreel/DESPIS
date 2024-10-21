values = [2.0725 2.0586 2.0341 2.002 1.9426 1.9202 1.9062 1.9045 1.9021 1.8960]*10^5;
x = [1 2 3 4 5 6 7 8 9 10];


figure; 
plot(x,values);
xlabel('distance (cm)');
ylabel('capacity (bit/sec)')
sgtitle('Signal capacity in function of distance from speaker')

 

 % Parameters
Fs = 16000;                % Sampling frequency (Hz)
f1 = 700;                  % Lower stopband frequency (Hz)
f2 = 3000;                 % Upper stopband frequency (Hz)
N = 100;                   % Filter order

% Normalize frequencies with respect to Nyquist frequency (Fs/2)
Wn = [f1 f2] / (Fs/2);

% Design band-stop filter using fir1
b = fir1(N, Wn, 'stop');

% Plot the frequency response
%freqz(b, 1, 1024, Fs);

% If you have a signal to filter, say x
y = filter(b, 1, wgn(Fs,1,0));    % Apply the filter to the signal
spectrogram(y,1024,512,1024,Fs);


