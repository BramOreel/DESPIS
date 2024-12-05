% Decision direction mode of adaptive linear filtering
%% Cleanup
clc; clear; close all;








%% Parameters.
M =4; % QAM constelllation size
Nq = log2(M);

nbQAMsymb =0; % Number of QAM symbols 
Hk = -0.7 + 0.8i; % Channel to consider (Only one frequency bin will be considered here)
alpha = 1; % Regularisation constant
SNR = 30; % Signal-to-noise-ratio [dB]

N_filter = 1; %number of tones in the filter





%% Construct data sequences.
nbBits = Nq*1000 ; % Number of bits required to generate nbQAMsymb symbols
bitstream = randi([0 1],nbBits,1); % bitstream of nbBits bits
Xk = qam_mod(bitstream,M); % QAM symbol sequence
Yk = Hk.*Xk; % REcorded QAM symbol sequence. This is our desired signal
Yk = awgn(Yk,SNR,"measured");

iN = 1; % Counter of stepsizes



%% Try to estimate wk as 1/Hk

for mu = [0.02 0.1 0.2 0.5] % List of stepsizes
    % NMLS filter implementation.
    % Initialise filters, reconstructed transmitted signal and error
    delta = 0.2;
    SNR = 30;
    w = zeros(length(Xk),1); 
    rec_Xk = zeros(1,length(Xk)); Ek = zeros(1,length(Xk));
    w(1) = (1 + delta)/conj(Hk); 
    for n = 1:length(Xk)-1
        % Apply filter.
        estXk = conj(w(n))*Yk(n);
        % Reconstruct transmitted signal.
        rec_Xk(n) = qammod(qamdemod(estXk,M),M);
        % Calculate error signal.
        Ek(n) =  estXk - rec_Xk(n) ;
        % Update filter.
        w(n+1) = w(n) + mu/(alpha + conj(Yk(n))*Yk(n))*Yk(n)*conj((Xk(n) - conj(w(n))*Yk(n)));
        if(abs(w(n+1)) > 10^30)
            break
        end
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

%We expect the value to ossicilate because our step size is rigid (should
%be dynamic)

%Think about question 5 i guess


end


