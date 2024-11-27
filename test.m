% Original frequency-domain signal (example, symmetric for real signals)
freq_domain = [1, 2, 3, 4, 0, 4, 3, 2]; % 8-point FFT with Hermitian symmetry
N = length(freq_domain); % Original number of points

% Step 1: Insert zeros between frequency samples
zero_padding_factor = 2; % Interpolation factor (e.g., 2x frequency resolution)
new_length = N * zero_padding_factor; % New length
freq_domain_interpolated = zeros(1, new_length); % Initialize zero-padded array

% Fill the first half of the spectrum with zeros in between
freq_domain_interpolated(1:zero_padding_factor:(N/2 + 1) * zero_padding_factor) = freq_domain(1:(N/2 + 1));

% Fill the second half (negative frequencies) of the spectrum, preserving Hermitian symmetry
freq_domain_interpolated(end - zero_padding_factor * (N/2 - 1):-zero_padding_factor:1 + zero_padding_factor) = ...
    freq_domain(end:-1:(N/2 + 2));

% Step 2: Inverse FFT to get the interpolated time-domain signal
time_domain_interpolated = ifft(freq_domain_interpolated, 'symmetric'); % Use 'symmetric' for real signals

% Optional: Check the interpolated frequency spectrum
freq_domain_check = fft(time_domain_interpolated);

% Visualization
disp('Original Frequency Domain:');
disp(freq_domain);

disp('Interpolated Frequency Domain:');
disp(freq_domain_interpolated);
