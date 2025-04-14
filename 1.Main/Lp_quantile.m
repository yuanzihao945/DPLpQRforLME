function xLp = Lp_quantile(x, p, taus, dec)
% Lp_quantile for Sample
%
% Input :
%   x       -  Sample Vector
%   p       -  the order of quantile, scalar
%   taus    -  Lp_quantile level
%   dec     -  Integer indicating the number of decimal places (round)
%
% Output:
%   xLp     -  Sample Lp_quantile
%
% Usage:
%   xLp = Lpquantile(x, p, taus, dec)
%           'taus' has default value '0 : 0.1 : 1';
%           'dec'  has default number '4';
%
% See also EVS, QUANTILE, QVS

% Author  : ZH.Yuan
% Update  : 2022/03/08 (First Version: 2022/03/03)
% Email   : zihaoyuan@whut.edu.cn (If any suggestions or questions)

if ~exist('taus', 'var') || isempty(taus)
    taus = 0 : 0.1 : 1;              % Setting default Lp_quantile level
end

if ~exist('dec', 'var') || isempty(dec)
    dec = 16;                       % Setting default decimal places
end

if ~isnumeric(x)
    error('observations are needed in vector form.')
end

if min(taus) < 0 || max(taus) > 1
    error('only asymmetries between 0 and 1 allowed.')
end

if p == 1

    xLp = quantile(x, taus);

elseif p == 2

    xLp = expectile(x, taus, dec);

else

    if min(size(x)) == 1

        xLpk = mean(x);
        xLp = 0 * taus;
        MaxTol = max(abs(x)) * 1e-06;

        for k = 1 : numel(taus)
            tau = taus(k);
            if tau == 0
                xLp(k) = min(x);
            elseif tau == 1
                xLp(k) = max(x);
            else
                for t = 1 : 200
                    d1 = sum(((x < xLpk) * (1 - tau) - (x >= xLpk) * tau) .* abs(x - xLpk).^(p - 1));
                    d2 = sum((p - 1) * abs(tau - (x - xLpk < 0)) .* abs(x - xLpk).^(p - 2));
                    de = d1 / d2;
                    xLpk = xLpk - de;
                    if max(abs(de)) < MaxTol
                        break
                    end
                end
                xLp(k) = xLpk;
            end
        end

        xLp = round(xLp, dec);

    else

        xLp = round(Lp_quantilemex(x, p, taus), dec);

    end

end


end



% function xLp = Lp_quantilepar(x, p, taus)
% P = size(x, 2);
% xLp = zeros(P, length(taus));
% parfor i = 1 : P
%     xi = x(:, i);
%     xLpk = mean(xi);
%     xLpv = 0 * taus;
%     MaxTol = max(abs(xi)) * 1e-06;
% 
%     for k = 1 : numel(taus)
%         tau = taus(k);
%         if tau == 0
%             xLpv(k) = min(xi);
%         elseif tau == 1
%             xLpv(k) = max(xi);
%         else
%             for t = 1 : 200
%                 d1 = sum(((xi < xLpk) * (1 - tau) - (xi >= xLpk) * tau) .* abs(xi - xLpk).^(p - 1));
%                 d2 = sum((p - 1) * abs(tau - (xi - xLpk < 0)) .* abs(xi - xLpk).^(p - 2));
%                 de = d1 / d2;
%                 xLpk = xLpk - de;
%                 if abs(de) < MaxTol
%                     break
%                 end
%             end
%             xLpv(k) = xLpk;
%         end
%     end
% 
%     xLp(i, :) = xLpv;
% 
% end
% end

