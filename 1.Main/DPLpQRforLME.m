function [beta, alpha, flambda_beta, flambda_alpha, lambda_beta_s, lambda_alpha_s, iter] = ...
    DPLpQRforLME(s, y, x, z, p, tau, opset)


%% Default optimize settings
if ~exist("opset", "var") || isempty(opset)
    opset = struct();
end

if ~exist("opset.eps", "var") || isempty(opset.eps)
    eps = 1e-4;
else
    eps = opset.eps;
end

if ~exist("opset.maxIter", "var") || isempty(opset.maxIter)
    maxIter = 1e2;
else
    maxIter = opset.maxIter;
end

if ~exist("opset.nfold", "var") || isempty(opset.nfold)
    nfold = 5;
else
    nfold = opset.nfold;
end

if ~exist("opset.criterion", "var") || isempty(opset.criterion)
    criterion = 'GACV';
else
    criterion = opset.criterion;
end

if ~exist("opset.nlambda", "var") || isempty(opset.nlambda)
    nlambda = 50;
else
    nlambda = opset.nlambda;
end

if ~exist("opset.threshbeta", "var") || isempty(opset.threshbeta)
    threshbeta = 1e-1;
else
    threshbeta = opset.threshbeta;
end

if ~exist("opset.threshalpha", "var") || isempty(opset.threshalpha)
    threshalpha = 1e-1;
else
    threshalpha = opset.threshalpha;
end


%% Dimension of data
n = length(unique(s));
N = length(y);
CVindex = CVgroup(N, nfold);

%% Extend the covariate associated with random effects to matrix
xe = [ones(N, 1) x];
ze = z;


%% Set the initial Value of alpha and beta
% Set alpha^{(0)} and calculate beta^{(0)} by solving the Lp_quantile Lasso regression
[beta_old, ~] = IterEstimation(xe, y, p, tau, nfold, nlambda, CVindex, criterion);
% Set initial difference of beta
beta_diff = 1;
% Set initial iter number
iter = 0;

%% The iterative algorithm can be restated for the Lasso
while beta_diff > eps

    iter = iter + 1;

    %% Update alpha
    % Calculate the modiﬁes residual r⌃{(m)}
    yres = y - xe * reshape(beta_old, length(beta_old), 1);
    % Calculate alpha^{(m+1)} by solving the Lp_quantile Lasso optimization
    [alpha, ~] = IterEstimationsupp(ze, yres, s, p, tau, nfold, nlambda, CVindex, criterion);

    %% Update beta
    % Calculate the new response ynew
    for is = 1 : n
        ynew(s == is, 1) = y(s == is) - ze(s == is, :) * alpha(is, :)';
    end
    % Calculate beta^{(m+1)} by solving the Lp_quantile Lasso optimization
    [beta_new, ~] = IterEstimation(xe, ynew, p, tau, nfold, nlambda, CVindex, criterion);

    %% Variable Selection
    alpha_select = mean(abs(alpha)) / sum(mean(abs(alpha)));
    index_alpha = alpha_select >= threshalpha;
    index_beta = abs(beta_new) >= threshbeta;
    % index_alpha = repmat((alpha_select / max(alpha_select)) >= threshalpha, n, 1);
    % index_beta = (abs(beta_new) / max(abs(beta_new))) >= threshbeta;
    beta_old_s = beta_new(index_beta);

    %% Update X and Z
    xenew = xe(:, index_beta);
    zenew = ze(:, index_alpha);
    
    %% Update alpha
    % Calculate the modiﬁes residual r⌃{(m)}
    yres = y - xenew * beta_old_s;
    % Calculate alpha^{(m+1)} by solving the quantile Lasso optimization
    [alpha_s, lambda_alpha_s(iter)] = IterEstimationsupp(zenew, yres, s, p, tau, nfold, nlambda, CVindex, criterion);

    %% Update beta
    % Calculate the new response ynew
    for is = 1 : n
        ynew(s == is, 1) = y(s == is) - zenew(s == is, :) * alpha_s(is, :)';
    end
    % Calculate beta^{(m+1)} by solving the quantile Lasso optimization
    [beta_new_s, lambda_beta_s(iter)] = IterEstimation(xenew, ynew, p, tau, nfold, nlambda, CVindex, criterion);

    %% Update beta and alpha
    beta_new(index_beta) = beta_new_s;
    beta_new(~index_beta) = 0;
    alpha(:, index_alpha) = alpha_s;
    alpha(:, ~index_alpha) = 0;

    %% Update beta_diff and iter number
    beta_diff = sum(abs(beta_new - beta_old)) / sum(beta_new ~= 0);
    beta_old = beta_new;
    
    %% Termination conditions
    if iter >= maxIter
        break
    end

end

beta = beta_old;
flambda_beta = lambda_beta_s(iter);
flambda_alpha = lambda_alpha_s(iter);

end









