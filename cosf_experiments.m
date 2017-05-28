function [regions, E_copula, u] = cosf_experiments(x,y)

u = pobs(x);
v = pobs(y);

[u,I] = sort(u);
v = v(I);

n = length(x); E_copula = zeros(1,n);
for ii=1:n
    E_copula(ii) = sum(u(ii)>=u & v(ii)>=v)/(n+1);
end

ss_minIdx = 1;
regions = {};
regionsIdx = 1;

for ii=2:n
    ss = E_copula(ss_minIdx:ii);
    if(~issorted(ss) && ~issorted(fliplr(ss)))      % just a way of testing if they are tied?
        
        regionVec = [ss_minIdx, ii];
        regions{regionsIdx} = regionVec; regionsIdx = regionsIdx + 1;
        ss_minIdx = ii;
        
    else
        if(ii==n)
            regionVec = [ss_minIdx, ii];
            regions{regionsIdx} = regionVec; regionsIdx = regionsIdx + 1;
        end
    end
end

end
