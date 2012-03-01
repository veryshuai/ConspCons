function [empty,out] = mpec_err(X,price,egr,v,w,scl,s,N,T,me_scale)

empty = []; %inequality constraints...mine are all equality
%tic 
warning off all;

param1 = X(1:29); %cobb douglas parameters
param2 = X(30:58); %stone geary parameters
alp = X(59); %weight of conspicuous consumption
vm = X(60); %v multiplier

for k = 1:N
v{k} = v{k}*vm; %scale v by the multiplier
end

upsiz = 29*2+2; %size of utility parameter vector
cfsiz = 28*3; %size of coefficient vector for each year-type

%set up param cells
cp = cell(T,N);
g_int = reshape(X(upsiz+1:end),N*cfsiz,T);
g = mat2cell(g_int,ones(1,N)*cfsiz,ones(1,T))';

options=optimset('Display','iter','jacobian','on','TolFun',1e-5,'TolX',1e-6,'TolCon',1e-6,'DerivativeCheck','off',...
            'GradObj','on','MaxIter',200,'LargeScale','off');

parfor t = 1:T
    for n = 1:N
        pscale = abs(util_FOC([param1,param2],g{t,n},w(:,t,n),price(t,:),scl,...
            g{t,n},alp,s{t,:},v{n},egr{t},1));
        cp{t,n} = fminunc(@(x) util_FOC([param1,param2],x,w(:,t,n),price(t,:),scl,...
            g{t,n},alp,s{t,:},v{n},egr{t},pscale),g{t,n},options);
    end
end

out_int = cell2mat(cp');
out = out_int(:) - X(upsiz+1:end);

out = out/me_scale;

%toc
%display(norm(out));

end