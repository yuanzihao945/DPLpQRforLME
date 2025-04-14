function value = criterionFun(criterion, N, SM, M)
switch criterion
    case 'SIC'
        value = log(SM ./ N) + M .* log(N) ./ (2 * N);
    case 'GACV'
        value = SM ./ (N - M);
end
value(isnan(value)) = inf;
end