function [theta, lambdaopt] = IterEstimation(x, y, p, tau, nfold, nlambda, CVindex, criterion)
N = length(y);
checkfun = @(x) abs(x).^p .* ((1 - tau) * (x < 0) + tau * (x >= 0));
SM_theta = zeros(nfold, nlambda);
M_theta = zeros(nfold, nlambda);
[~, fit] = lasso(x, y, 'NumLambda', nlambda, 'Intercept', false, 'Standardize', false);
lambdas = fit.Lambda;

for infold = 1 : nfold
    xi = x(CVindex ~= infold, :);
    yi = y(CVindex ~= infold);
    xri = x(CVindex == infold, :);
    yri = y(CVindex == infold);
    if p == 1
        theta_i = lasso(xi, yi, 'Lambda', lambdas, 'Intercept', false, 'Standardize', false);
    else
        theta0 = repmat(pinv(xi' * xi) * xi' * yi, 1, nlambda);
        theta_i = Theta_Fastmex(yi, xi, lambdas, theta0, p, tau, 1e-4);
    end
    SM_theta(infold, :) = sum(checkfun(yri - xri * theta_i));
    M_theta(infold, :) = sum(theta_i ~= 0);
end

Cvalue_theta = sum(criterionFun(criterion, N, SM_theta, M_theta));
Cvalue_theta(isnan(Cvalue_theta)) = inf;
optp = find(Cvalue_theta == min(Cvalue_theta));
lambdaopt = lambdas(optp(1));
if p == 1
    theta = theta_i(:, optp(1));
else
    theta = Theta_Fastmex(y, x, lambdaopt, theta_i(:, optp(1)), p, tau, 1e-4);
end

end