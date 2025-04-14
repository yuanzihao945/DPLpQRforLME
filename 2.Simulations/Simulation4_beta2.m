%% Clean the workspace
clearvars
close all
clc

%% Set Parameters
n = 10;                     % Object Number
m = 10;                     % Time point length                                             
N = n * m;                  % Sample Number
beta0 = [0; 3; 4; 0; 0; 1.5; 2; ones(2993, 1)];   % Coefficient of fixed effect
k = size(beta0, 1);          % Number of the fixed effect variable
rho = 0.5;                  % Correlation coefficient of fixed effect variables
sigmaX = rho.^(abs(((1:k) - (1:k)')));
s = repmat((1:n)', m, 1);   % Object Indicator variable 

P = 6;                      % Number of the Random effect variable
sigmaZ = rho.^(abs(((1:P) - (1:P)')));
sigmaalpha = diag([ones(1, P/2),  zeros(1, P/2)]);   % Covariance matrix of random effects
sigma = 0.5;

for tau = [0.25 0.5 0.75]
    for ip = 1 : 11
        p = 1 + (ip - 1) / 10;
        for irepeat = 1 : 100
            irun = 1;
            while irun
                x = mvnrnd(zeros(1, k), sigmaX, N);          % Generation of fixed effect variables
                z = mvnrnd(zeros(1, P), sigmaZ, N);          % Generation of Random effect variables
                alpha0 = mvnrnd(zeros(1, P), sigmaalpha, n); % Coefficient of Random effect
                y = zeros(N, 1);
                for i = 1 : N
                    y(i) = 1.5 + x(i, :) * beta0 + z(i, :) * alpha0(s(i), :)' + sigma * randn(1);
                end
                try
                    [beta, alpha, lambda_beta, lambda_alpha, lambda_beta_s, lambda_alpha_s] ...
                        = DPLpQRforLME(s, y, x, z, p, tau);
                    irun = 0;
                catch
                    irun = 1;
                end
            end
            betahat(:, irepeat) = beta;
            MSE(irepeat) = (beta(2:end) - beta0)' * pinv(sigmaX) * (beta(2:end) - beta0);
            disp(['tau_', num2str(tau), ...
                '_L_', num2str(p), '_', num2str(irepeat), '-th done.'])
        end

        MSEm(ip, :) = mean(MSE);
        MSEs(ip, :) = std(betahat,[], 2);
        Corr = sum(((betahat ~= 0) + [1; beta0 ~= 0]) == 2, 1);
        Incorr = sum(((betahat ~= 0) + [0; beta0 == 0]) == 2, 1);
        Corrm(ip, :) = mean(Corr);
        Corrs(ip, :) = std(Corr);
        Incorrm(ip, :) = mean(Incorr);
        Incorrs(ip, :) = std(Incorr);
        CorrNum(ip, :) = sum(((betahat ~= 0) + [1; beta0 ~= 0]) == 2, 2)' ...
            + sum(((betahat == 0) + [0; beta0 == 0]) == 2, 2)';
    end
    writematrix([MSEm, MSEs, Corrm, Corrs, Incorrm, Incorrs, CorrNum], ...
        ['Simulation4_n10m10_beta2_tau', num2str(tau), '.csv'])
end







