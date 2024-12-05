% Decision direction mode of adaptive linear filtering
%% Cleanup
clc; clear; close all;








%% Parameters.
M =8; % QAM constelllation size
Nq = log2(M);
nbQAMsymb =0; % Number of QAM symbols 
Hk = -1.3 + 5i; % Channel to consider (Only one frequency bin will be considered here)
alpha = 0; % Regularisation constant
SNR = 30; % Signal-to-noise-ratio [dB]

%% Construct data sequences.
nbBits =0 ; % Number of bits required to generate nbQAMsymb symbols
bitstream = randi([0 1],Nq*1000,1); % bitstream of nbBits bits

%This corresponds to the transmitted QAM symbol sequence Xk at a single tone k.
Xk = qam_mod(bitstream,M); % QAM symbol sequence 


Yk = Hk.*Xk; % REcorded QAM symbol sequence
iN = 1; % Counter of stepsizes

Yk = awgn(Yk,SNR,"measured");



%% biem

for mu = [0.5; 1] % List of stepsizes
    % NMLS filter implementation.
    % Initialise filters, reconstructed transmitted signal and error
    
    w = zeros(length(Xk),1); 
    rec_Xk = zeros(1,length(Xk)); Ek = zeros(1,length(Xk));
    delta = 0.2;
    w(1) = (1+delta)/conj(Hk); 
    for n = 1:length(Xk)-1
        % Apply filter.
        estXk = ;
        % Reconstruct transmitted signal.
        rec_Xk(n) = ;
        % Calculate error signal.
        Ek(n) = ;
        % Update filter.
        w(n+1) = w(n)+mu/(alpha + );
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

