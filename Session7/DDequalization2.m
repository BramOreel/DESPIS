% Decision direction mode of adaptive linear filtering
%% Cleanup
clc; clear; close all;

%% Parameters.
N_q = 3;
M = 2^N_q; % QAM constelllation size
nbQAMsymb = 1000; % Number of QAM symbols 
Hk = complex(0.541, 0.021); % Channel to consider (Only one frequency bin will be considered here)
delta = 0.2; 
alpha = 1; % Regularisation constant
SNR = 30; % Signal-to-noise-ratio [dB]

%% Construct data sequences.
nbBits = nbQAMsymb*N_q; % Number of bits required to generate nbQAMsymb symbols
bitstream = randi([0,1], [nbBits, 1]); % bitstream of nbBits bits
Xk = qam_mod(bitstream,M); % QAM symbol sequence
Yk = Xk*Hk; % REcorded QAM symbol sequence
Yk = awgn(Yk,SNR,"measured");
iN = 1; % Counter of stepsizes

for mu = [0.02, 0.1, 0.2, 0.5, 1] % List of stepsizes
    % NMLS filter implementation.
    % Initialise filters, reconstructed transmitted signal and error
    w = zeros(1,nbQAMsymb); 
    rec_Xk = zeros(1,nbQAMsymb); Ek = zeros(1,nbQAMsymb);
    w(1) = (1 + delta)*1/conj(Hk); 
    for n = 1:length(Xk)
        % Apply filter.
        estXk = Yk(n)*conj(w(n));
        % Reconstruct transmitted signal.
        rec_Xk(n) = qammod(qamdemod(estXk, M, UnitAveragePower=true),M,UnitAveragePower=true,PlotConstellation=false);
        % Calculate error signal.
        Ek(n) = rec_Xk(n) - estXk ;
        % Update filter.
        w(n + 1) = w(n) + mu * (Yk(n) / (alpha + conj(Yk(n)) * Yk(n))) * conj(Ek(n)); %d(L+1): desired signal
    end

    %BER = ber(rec_Xk, Xk)
    legendCell{iN} = num2str(mu,'mu=% .2f');
    iN = iN + 1;

% Plot results.
figure(1);
semilogy(abs(conj(w) - 1/Hk));
title('NMLS estimation error.');
xlabel('Iteration'); ylabel('Error magnitude');
hold on
legend(legendCell)


end


