function xev = expectile(x, tau, dec)
% Expectile for Sample
% 
% Input :
%   x       -  Sample Vector
%   tau     -  Expectile level
%   dec     -  Integer indicating the number of decimal places (round) 
%
% Output:  
%   xev     -  Sample Expectile
%
% Usage:
%   xev = EXPECTILE(x, tau, dec)
%           'tau' has default value '0 : 0.1 : 1';
%           'dec'  has default number '4';
%
% See also EVS, QUANTILE, QVS

% Author  : ZH.Yuan
% Update  : 2022/03/08 (First Version: 2022/03/03)
% Email   : zihaoyuan@whut.edu.cn (If any suggestions or questions)

if ~exist('tau', 'var') || isempty(tau)
    tau = 0 : 0.1 : 1;              % Setting default expectile level
end

if ~exist('dec', 'var') || isempty(dec)
    dec = 16;                       % Setting default decimal places
end

if ~isnumeric(x)
    error('observations are needed in vector form.')
end

if min(tau) < 0 || max(tau) > 1
    error('only asymmetries between 0 and 1 allowed.')
end

[n, p] = size(x);

if min(n, p) == 1

    xek = mean(x);
    xev = 0 * tau;
    MaxTol = max(abs(x)) * 1e-06;

    for k = 1 : numel(tau)
        p = tau(k);
        if p == 0
            xev(k) = min(x);
        elseif p == 1
            xev(k) = max(x);
        else
            for t = 1 : 200
                w = abs(p - (x < xek));
                xek_new = sum(w .* x) / sum(w);
                dxek = max(abs(xek_new - xek));
                xek = xek_new;
                if dxek < MaxTol
                    break
                end
            end
            xev(k) = xek;
        end
    end

    xev = round(xev, dec);

else

    xev = round(expectilemex(x, tau), dec);

end

end


% %%
% function xevm = expectilepar(x, tau)
% 
% P = size(x, 2);
% xevm = zeros(P, length(tau));
% parfor iP = 1 : P
%     xiP = x(:, iP);
%     xek = mean(xiP);
%     xev = 0 * tau;
%     MaxTol = max(abs(xiP)) * 1e-06;
% 
%     for k = 1 : numel(tau)
%         p = tau(k);
%         if p == 0
%             xev(k) = min(xiP);
%         elseif p == 1
%             xev(k) = max(xiP);
%         else
%             for t = 1 : 200
%                 w = abs(p - (xiP < xek));
%                 xek_new = sum(w .* xiP) / sum(w);
%                 dxek = max(abs(xek_new - xek));
%                 xek = xek_new;
%                 if dxek < MaxTol
%                     break
%                 end
%             end
%             xev(k) = xek;
%         end
%     end
% 
%     xevm(iP, :) = xev;
% 
% end
% end
% 
% 
