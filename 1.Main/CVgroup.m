%% Function for CV
function CVlist = CVgroup(datasize, nfold)
if ~exist("nfold", "var") || isempty(nfold)
    nfold = 5;
end
splace = ceil(datasize * (1 : nfold) / nfold);
label = randperm(datasize);
CVlist(label(1 : splace(1))) = 1;
for i = 2 : nfold
    CVlist(label((splace(i-1)+1) : splace(i))) = i;
end
end