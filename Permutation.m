function [tw, perm, bags] = Permutation(nodeL,fnodeL,tnodeL,alpha)

    nb = length(nodeL);
    nl = length(fnodeL);
    
    busN = sparse(nodeL, ones(nb,1), (1:nb) );
    fbusN = busN(fnodeL);
    tbusN = busN(tnodeL);    

    perm = zeros(nb,1);
    neighbour = cell(nb,1);
    bags = cell(nb,1);

    for ii = 1 : nl
        neighbour{fbusN(ii)} = union( neighbour{fbusN(ii)}, tnodeL(ii) );
        neighbour{tbusN(ii)} = union( neighbour{tbusN(ii)}, fnodeL(ii) );
    end

    fil = Fill_In(busN,neighbour,nodeL);
    deg = Deg(busN,neighbour,nodeL);

    tw = 1;
    for ii = 1 : nb
%         ii
        [a,m] = min(fil);
        if a ~= 0
            [a,m] = min( alpha*deg + fil);
        end

        perm(ii) = nodeL(m);
        bb = union(perm(ii),neighbour{m});

        bags{ii} = bb;
        if  length(bb) > tw + 1
            tw = length(bb) - 1;
        end

        nn = busN(neighbour{m});
        nn2 = neighbour{m};
        for jj = 1 : length(nn)
            kk = nn(jj);
            nn2 = union(nn2,neighbour{kk});
            neighbour{kk} = setdiff(neighbour{kk},nodeL(m));
            neighbour{kk} = union(neighbour{kk},setdiff(neighbour{m},nodeL(kk)));
        end

        ff = Fill_In(busN,neighbour,nn2);
        fil(busN(nn2)) = ff;

        dd = Deg(busN,neighbour,neighbour{m});
        deg(nn) = dd;    

        fil(m) = Inf;
        deg(m) = Inf;

        neighbour{m} = 0;
    end
%     disp('Upper bound on treewidth:');
%     disp(tw);
end