% Parameters
T = 1;  % Time window duration
f_max = 5;  % Maximum frequency for plotting
fs = 1000;  % Sampling frequency for better resolution

% Time domain - Rectangular window
t = -T:1/fs:T;  % Time vector for one period of the rectangular window
rect_window = double(abs(t) <= T/2);  % Rectangular window function

% Frequency domain - Sinc of the rectangular window (Fourier transform)
f = -f_max:1/fs:f_max;  % Frequency vector for plotting
rect_sinc = T * sinc(T * f);  % Fourier transform of rectangular window

% Time domain - Sinc function (convolution with rectangular window in freq domain)
sinc_t = sinc(t);  % Sinc function

% Frequency domain - Convolution result
conv_result = conv(rect_sinc, sinc_t, 'same');  % Convolution in the frequency domain

% Plotting
figure;
subplot(2,2,1);
plot(t, rect_window);
title('Rectangular Window in Time Domain');
xlabel('Time');
ylabel('Amplitude');
axis([-T T -0.2 1.2]);

subplot(2,2,2);
plot(f, rect_sinc);
title('Sinc Function in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(2,2,3);
plot(t, sinc_t);
title('Sinc Function in Time Domain');
xlabel('Time');
ylabel('Amplitude');

subplot(2,2,4);
plot(f, conv_result);
title('Convolution of Rectangular Window and Sinc in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
axis([-f_max f_max -1 1]);

