function r = resid(mc,ac,egc)
%This function takes model consumption, actual consumption, and actual
%expenditure and expenditure grid cutoffs and returns a likelihood

%make a big model consumption grid (note that we leave out the really rich
%guys...we don't consider them below anyway)
bmcg = cell(size(ac));
r = 0;
obs_num = 0;
sq_r = 0;
for j=1:4
    for m = 1:18
        bmcg{m,j} = zeros(size(ac{m,j},1),29);
        obin = 1;
        for k = 1:size(mc,1)
            bmcg{m,j}(obin:egc{m,j}(k),:) = repmat(mc(k,:,m,j),egc{m,j}(k)-obin+1,1);
            obin = egc{m,j}(k);
        end
        bmcg{m,j} = bmcg{m,j}(egc{m,j}(1):egc{m,j}(end),:); %shrink out extremes
        ac{m,j} = ac{m,j}(egc{m,j}(1):egc{m,j}(end),:); %shrink out extremes
        obs_num = obs_num + size(ac{m,j},1) * size(ac{m,j},2);
        sq_r = sum(sum((ac{m,j}-bmcg{m,j}).^2));
    end
end

var_est = sq_r/obs_num;
display(var_est);

for j=1:4
    for m = 1:18
        r = r + sum(sum(log(normpdf(ac{m,j},bmcg{m,j},var_est*ones(size(ac{m,j})))))); %recursively calculate likelihood
    end
end

end