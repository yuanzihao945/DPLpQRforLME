function [theta, lambdaopt] = IterEstimationsupp(x, y, s, p, tau, nfold, nlambda, CVindex, criterion)
N = length(y);
checkfun = @(x) abs(x).^p .* ((1 - tau) * (x < 0) + tau * (x >= 0));
SM_theta = zeros(nfold, nlambda);
M_theta = zeros(nfold, nlambda);


[~, fit] = lasso(x, y, 'NumLambda', nlambda, 'Intercept', false, 'Standardize', false);
lambdas = fit.Lambda;

for infold = 1 : nfold
    for is = 1 : length(unique(s))
        xi = x((CVindex ~= infold)' & (s == is), :);
        yi = y((CVindex ~= infold)' & (s == is));
        xri = x((CVindex == infold)' & (s == is), :);
        yri = y((CVindex == infold)' & (s == is));
        if p == 1
            theta_i = lasso(xi, yi, 'Lambda', lambdas, 'Intercept', false, 'Standardize', false);
        else
            theta0 = repmat(pinv(xi' * xi) * xi' * yi, 1, nlambda);
            theta_i = Theta_Fastmex(yi, xi, lambdas, theta0, p, tau, 1e-4);
        end
        SM_theta(infold, :) = SM_theta(infold, :) + sum(checkfun(yri - xri * theta_i));
        M_theta(infold, :) = M_theta(infold, :) + sum(theta_i ~= 0);
    end
end

Cvalue_theta = sum(criterionFun(criterion, N, SM_theta, M_theta / nfold));
Cvalue_theta(isnan(Cvalue_theta)) = inf;
optp = find(Cvalue_theta == min(Cvalue_theta));
lambdaopt = lambdas(optp(1));
for is = 1 : length(unique(s))
    xs = x(s == is, :);
    ys = y(s == is);
    if p == 1
        theta(is, :) = lasso(xs, ys, 'Lambda', lambdaopt, 'Intercept', false, 'Standardize', false);
    else
        theta0 = pinv(xs' * xs) * xs' * ys;
        theta(is, :) = Theta_Fastmex(ys, xs, lambdaopt, theta0, p, tau, 1e-4);
    end
end

end