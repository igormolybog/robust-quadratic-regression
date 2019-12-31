function [V] = m2v(M,ndx)
    n = length(ndx);
    MM = M - M.*sparse(1:n,1:n,1-1/sqrt(2),n,n);
    V = sqrt(2) * [real(MM(triu(ndx)~=0)); imag(MM(triu(ndx,1)~=0))];
end
