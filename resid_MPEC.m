function r = resid_MPEC(x,ac,ce,scl,N,T,resid_scale)
%This function takes model coefficients, actual budget shares and total
%expenditure level

upsiz = 29*2+2; %size of utility parameter vector
cfsiz = 28*3; %size of coefficient vector for each year-type
obs_num = 0;
sq_r = 0;
bmcg = cell(size(ac)); %big model consumption grid
tic
for t = 1:T
    for n = 1:N
        cf = x(upsiz+cfsiz*(N*(t-1)+(n-1))+1:upsiz+cfsiz*(N*(t-1)+n),1);
        cf = reshape(cf,size(cf,1)/28,28)';
        for k = 1:size(cf,2)
            cf(:,k) = cf(:,k)/scl^(k-1);
        end

        %get budget shares from coefficients        
        bmcg{t,n} = bs(cf,ce{t,n});
        sq_r = sq_r + sum(sum((ac{t,n}-bmcg{t,n}).^2));
        obs_num = obs_num +size(ac{t,n},1);
    end
end
%toc
var_est = sq_r/obs_num;
%display(var_est);

r = 0;
for n=1:N
    for t = 1:T
        r = r + -sum(sum(log(normpdf(ac{t,n},bmcg{t,n},var_est*ones(size(ac{t,n})))))); %recursively calculate likelihood
    end
end

r = r/resid_scale;

% fid = fopen('2-18-2012-fp.txt','a');
% fprintf(fid,'resid\n\n');
% fprintf(fid,'%f\n',r);

end