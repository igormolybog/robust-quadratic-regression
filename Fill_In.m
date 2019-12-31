function [ff] = Fill_In(busN,neighbour,bb)
    ff = zeros(length(bb),1);
    for ii = 1 : length(bb)
        qq = busN(bb(ii));
        dd = length(neighbour{qq});
        nn = busN(neighbour{qq});
        ll = 0;
        for jj = 1 : dd
            kk = nn(jj);
            ll = ll + length(intersect(neighbour{qq}, neighbour{kk}));
        end
        ff(ii) = (dd*(dd - 1) - ll)/2;
    end
end