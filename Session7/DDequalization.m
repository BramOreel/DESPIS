% Decision direction mode of adaptive linear filtering
%% Cleanup
clc; clear; close all;

%% Parameters.
M =; % QAM constelllation size
nbQAMsymb =; % Number of QAM symbols 
Hk = ; % Channel to consider (Only one frequency bin will be considered here)
alpha = ; % Regularisation constant
SNR = ; % Signal-to-noise-ratio [dB]

%% Construct data sequences.
nbBits = ; % Number of bits required to generate nbQAMsymb symbols
bitstream = randi(); % bitstream of nbBits bits
Xk = qam_mod(); % QAM symbol sequence
Yk = ; % REcorded QAM symbol sequence
iN = 1; % Counter of stepsizes

for mu = [] % List of stepsizes
    % NMLS filter implementation.
    % Initialise filters, reconstructed transmitted signal and error
    w = zeros(); 
    rec_Xk = zeros(); Ek = zeros();
    w(1) = ; 
    for n = 1:length(Xk)-1
        % Apply filter.
        estXk = ;
        % Reconstruct transmitted signal.
        rec_Xk(n) = ;
        % Calculate error signal.
        Ek(n) = ;
        % Update filter.
        w(n+1) = ;
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


