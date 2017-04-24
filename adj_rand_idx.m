function [ARI] = adj_rand_idx(naiItr,commonSet)

% make contingency
nPartitions = length(naiItr);
ARI = nan(nPartitions);
for a = 1:nPartitions-1
    for b = a+1:nPartitions
        p1 = naiItr{a};
        p2 = naiItr{b};
        lenI = length(p1);
        lenJ = length(p2);
        conTable = nan(lenI,lenJ);
        binomTable = nan(lenI,lenJ);
        for I = 1:lenI
            for J = 1:lenJ
                conTable(I,J) = length(intersect(p1{I},p2{J}));
                if conTable(I,J)>2
                    binomTable(I,J) = nchoosek(conTable(I,J),2);
                else
                    binomTable(I,J) = conTable(I,J);
                end
            end
        end
        aBinom = nan(1,lenI);
        for I = 1:lenI
            iSum = sum(conTable(I,:));
            if iSum>1
                aBinom(1,I) = nchoosek(iSum,2);
            else
            	aBinom(1,I) = iSum;
            end
        end
        
        bBinom = nan(lenJ,1);
        for J = 1:lenJ
            jSum = sum(conTable(:,J));
            if jSum>1
                bBinom(J,1) = nchoosek(jSum,2);
            else
                bBinom(J,1) = jSum;
            end
        end
        
        aSums = sum(aBinom);
        bSums = sum(bBinom);
        n = length(commonSet);
        nBinom = nchoosek(n,2);
        ARInumerator = sum(sum(binomTable))-((aSums*bSums)/nBinom);
        ARIdenominator = .5*(aSums+bSums)-((aSums*bSums)/nBinom);
        
        ARI(a,b) = ARInumerator/ARIdenominator;
        ARI(b,a) = ARI(a,b);
    end
end
1;