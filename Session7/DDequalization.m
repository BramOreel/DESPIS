% Decision direction mode of adaptive linear filtering
%% Cleanup
clc; clear; close all;








%% Parameters.
M =16; % QAM constelllation size
Nq = log2(M);
constellation = (-Nq+1:2:Nq-1) + 1j*(-Nq+1:2:Nq-1)';
nbQAMsymb =0; % Number of QAM symbols 
Hk = -1.3 + 5i; % Channel to consider (Only one frequency bin will be considered here)
alpha = 1; % Regularisation constant
SNR = 30; % Signal-to-noise-ratio [dB]

N_filter = 1; %number of tones in the filter


%% Construct data sequences.
nbBits = Nq*1000 ; % Number of bits required to generate nbQAMsymb symbols
bitstream = randi([0 1],nbBits,1); % bitstream of nbBits bits
Xk = qam_mod(bitstream,M); % QAM symbol sequence
Yk = Hk.*Xk; % REcorded QAM symbol sequence. This is our desired signal
Yk_noise = awgn(Yk,SNR,"measured");

iN = 1; % Counter of stepsizes
legendCell = {};



%% Try to estimate wk as 1/Hk

for mu = [0.02, 0.1, 0.2, 0.5, 1, 5, 10] % List of stepsizes
    % NMLS filter implementation.
    % Initialise filters, reconstructed transmitted signal and error
    delta = 0.2;
    SNR = 30;
    w = zeros(length(Xk),1); 
    rec_Xk = zeros(1,length(Xk)); Ek = zeros(1,length(Xk));
    w(1) = (1 + delta)/conj(Hk); 
    for n = 1:length(Xk)-1
        % Apply filter.
        estXk = w(n)'*Yk_noise(n);
        % Reconstruct transmitted signal.
        rec_Xk(n) = qammod(qamdemod(estXk,M),M);
        % Calculate error signal.
        Ek(n) =  Xk(n) - rec_Xk(n) ;
        % Update filter.
        %w(n+1) = w(n) + mu/(alpha + conj(Yk_noise(n))*Yk_noise(n))*Yk_noise(n)*conj((Xk(n) - conj(w(n))*Yk_noise(n)));
        w(n+1) = w(n) + mu/(alpha + conj(Yk_noise(n+1))*Yk_noise(n+1))*Yk_noise(n+1)*conj((Xk(n+1) - conj(w(n))*Yk_noise(n+1)));
    end
    legendCell{iN} = num2str(mu,'mu=% .2f');
    iN = iN + 1;

% Plot results.
figure(1);
semilogy(abs(1./w'-Hk));
title('NMLS estimation error.');
xlabel('Iteration'); ylabel('Error magnitude');
hold on
legend(legendCell)


end


