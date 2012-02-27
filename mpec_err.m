function [empty,out] = mpec_err(X,price,egr,v,w,scl,s,N,T)

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

g = cell(T,N);
cp = cell(T,N);
%set up param cells
for t = 1:T
    for n = 1:N
        g{t,n} = X(upsiz+cfsiz*(N*(t-1)+(n-1))+1:upsiz+cfsiz*(N*(t-1)+n),1);
    end
end

options=optimset('Display','off','jacobian','on','TolFun',1e-5,'TolX',1e-5,'TolCon',1e-6,'DerivativeCheck','off',...
            'GradObj','on','MaxIter',200,'LargeScale','off','Algorithm','active-set');
        
parfor t = 1:T
    for n = 1:N
        cp{t,n} = fminunc(@(x) util_FOC([param1,param2],x,w(:,t,n),price(t,:),scl,...
            g{t,n},alp,s{t},v{n},egr{t}),g{t,n},options);
    end
end

out = zeros(T*N*28*3,1);
for t = 1:T
    for n = 1:N
        out(cfsiz*(N*(t-1)+(n-1))+1:cfsiz*(N*(t-1)+n),1) = cp{t,n}-g{t,n};
    end
end

%toc
%display(norm(out));

end