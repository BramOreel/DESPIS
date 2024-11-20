function [H_est] = ofdm_channelest_pilots(received_pilots, transmitted_pilots, N, L)
% Estimates the channel frequency response using pilot tones.
%
% INPUT:
% received_pilots     Mx1     Received pilot symbols in the frequency domain.
% transmitted_pilots  Mx1     Transmitted pilot symbols in the frequency domain.
% N                   1x1     Number of OFDM tones (FFT size).
% L                   1x1     Length of the channel impulse response (equal to CP length).
%
% OUTPUT:
% H_est               Nx1     Estimated channel frequency response.

%% Define pilot positions
% Use odd tones (1, 3, 5, ..., N-1) as pilot tones.
pilot_indices = 1:2:N; % Odd indices for pilot tones

% Calculate the channel frequency response at pilot positions
H_pilots = received_pilots ./ transmitted_pilots; % Element-wise division

%% Interpolation to estimate the full channel frequency response

H_est = zeros(N, 1);

H_est(pilot_indices) = H_pilots;

h_time = ifft(H_est);

h_truncated = h_time;
h_truncated(L+1:end) = 0;

H_est = fft(h_truncated);

end