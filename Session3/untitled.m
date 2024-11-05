% Parameters
numSubcarriers = 8;       % Small number of subcarriers to observe overlap
M = 16;                   % Modulation order (e.g., 16-QAM)
symbolRate = 1e3;         % Symbol rate (1 kHz symbol duration)
fs = numSubcarriers * symbolRate;  % Sampling frequency

% Generate random QAM symbols for each subcarrier
dataSymbols = randi([0 M-1], numSubcarriers, 1);
qamSymbols = qammod(dataSymbols, M, 'UnitAveragePower', true);

% Apply IFFT to convert from frequency to time domain
timeDomainSignal = ifft(qamSymbols);

% Plotting the frequency domain and time-domain signals
figure;
subplot(2,1,1);
stem(0:numSubcarriers-1, abs(qamSymbols)); 
title('Frequency-Domain (QAM symbols on subcarriers)');
xlabel('Subcarrier Index');
ylabel('Magnitude');

subplot(2,1,2);
plot(real(timeDomainSignal)); 
title('Time-Domain Signal after IFFT (Overlapping Carriers)');
xlabel('Sample Index');
ylabel('Amplitude');
