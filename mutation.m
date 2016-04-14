function pnew = mutation(p,dims,currentgen,totgen,mutrate,lowerbd,upperbd)
% mutation operation
% By: Free Xiong; 2014-12-10
x = (1-currentgen/totgen)^(5/mutrate);
if x>1e-4
    whichdim = randperm(dims,1);
    mutrange = (upperbd(whichdim)-lowerbd(whichdim))*x;
    ub = p(whichdim)+mutrange;
    lb = p(whichdim)-mutrange;
    p(whichdim) = lb+rand()*(ub-lb);
end
pnew = p;