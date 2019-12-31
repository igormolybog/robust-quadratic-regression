function [M] = v2m(V,ndx)    
    n = length(ndx);
    l = length(V);
    d = trace(ndx);
    M = sparse(n,n);
    
    [fxd,txd] = find(triu(ndx));
    [fx,tx] = find(triu(ndx,1));
    
%     M = M + sparse(fxd, txd, V(1:(l+n)/2), n, n)/sqrt(2);
%     M = M + sparse(txd, fxd, V(1:(l+n)/2), n, n)/sqrt(2);
%     M = M + 1i*sparse(fx, tx, V(1+(l+n)/2:l), n, n)/sqrt(2);
%     M = M - 1i*sparse(tx, fx, V(1+(l+n)/2:l), n, n)/sqrt(2);
    M = M + sparse(fxd, txd, V(1:(l+d)/2), n, n)/sqrt(2);
    M = M + sparse(txd, fxd, V(1:(l+d)/2), n, n)/sqrt(2);
    M = M + 1i*sparse(fx, tx, V(1+(l+d)/2:l), n, n)/sqrt(2);
    M = M - 1i*sparse(tx, fx, V(1+(l+d)/2:l), n, n)/sqrt(2);
    M = M + M .* sparse(1:n,1:n,1/sqrt(2)-1,n,n);
end
