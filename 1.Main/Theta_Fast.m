function thetahat = Theta_Fast(yi, xi, lambdas, theta0, p, tau, eps)

nlambda = numel(lambdas);
thetahat = zeros(size(theta0, 1), nlambda);
parfor ilambda = 1 : nlambda
    iter = 1;
    theta = theta0(:, ilambda);
    lambda = lambdas(ilambda);
    lambdap = 1;
    
    while iter < 200
        delta = yi - xi * theta;
        checkfun = abs(delta).^p .* ((1 - tau) * (delta < 0) + tau * (delta >= 0));
        FX = sum(checkfun) + lambda * sum(abs(theta));
        dcheckfun = p * sign(delta) .* abs(delta).^(p-1) ...
            .* ((1 - tau) * (delta < 0) + tau * (delta >= 0));
        DFX = - xi' * dcheckfun + lambda * sum(sign(theta));
        d2checkfun = p * (p-1) .* abs(delta).^(p-2) ...
            .* ((1 - tau) * (delta < 0) + tau * (delta >= 0));
        D2FX = (xi .* d2checkfun)' * xi;
        FXnew = FX + 1;

        while FXnew > FX
            epsilon = - pinv(D2FX + lambdap) * DFX;
            thetanew = theta - epsilon;
            deltanew = yi - xi * thetanew;
            checkfunnew = abs(deltanew).^p .* ((1 - tau) * (deltanew < 0) + tau * (deltanew >= 0));
            FXnew = sum(checkfunnew) + lambda * sum(abs(thetanew));
            if FXnew > FX
                lambdap = lambdap / 0.1;
            else
                lambdap = lambdap * 0.1;
                theta = thetanew;
            end
        end
        
        if (abs(FX - FXnew) / abs(FX)) < eps
            break
        end
        iter = iter + 1;
        
    end
    thetahat(:, ilambda) = theta;
end

end
