function [ff] = Deg(busN,neighbour,bb)
    ff = zeros(length(bb),1);
    for ii = 1 : length(bb)
        qq = busN(bb(ii));
        ff(ii) = length(neighbour{qq});
    end
end