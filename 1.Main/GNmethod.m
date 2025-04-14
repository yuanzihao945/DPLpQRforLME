function Xopt = GNmethod(X0, DX, D2X)
iter = 1;
while iter < 200
    diff = pinv(D2X(X0)) * DX(X0);
    X = X0 - diff;
    if max(sum(abs(diff)) / sum(abs(X0))) < 1e-3
        break
    end
    iter = iter + 1;
    X0 = X;
end
Xopt = X;
end