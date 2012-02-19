function [empty,dis] = mpec_err(X,price,egr,v,w,scl)

empty = []; %inequality constraints...mine are all equality

warning off all;

%tic

%diary off
%diary 2-4-12.txt


    
param1 = X(1:29); %cobb douglas parameters
param2 = X(30:58); %stone geary parameters
alp = X(59); %weight of conspicuous consumption
vm = X(60); %v multiplier

%alp = 0;

v = v*vm; %scale v by the multiplier

upsiz = 29*2+2; %size of utility parameter vector
cfsiz = 28*3; %size of coefficient vector for each year-type

dis = zeros(28*3*18*4,1); %this holds the distance between the guessed coefficients and optimized coefficients
for t = 1:18
    for n = 1:4
        %display(4*(t-1)+(n-1)+1);
        %sobol points
        s = sobolset(1);
        s = net(s,10);
        s = egr(t,1) + (egr(t,end)-egr(t,1))*s;
        s = sortrows(s);

        options=optimset('Display','off','jacobian','on','TolFun',5e-3,'TolX',5e-3,'TolCon',1e-6,'DerivativeCheck','off',...
            'GradObj','on','MaxIter',100,'LargeScale','off');%,'Algorithm','active-set');
        
        g = X(upsiz+cfsiz*(4*(t-1)+(n-1))+1:upsiz+cfsiz*(4*(t-1)+n),1);
        
        cp = fminunc(@(x) util_FOC([param1,param2],x,w(:,t,n),price(t,:),scl,...
            g,alp,s,v(n,:),egr(t,:)),g,options);

        dis(cfsiz*(4*(t-1)+(n-1))+1:cfsiz*(4*(t-1)+n),1) = cp-g;
    end
end
%toc

fid = fopen('2-18-2012-fp.txt','a');
fprintf(fid,'eq_err\n\n');
fprintf(fid,'%f\n',norm(dis)/norm(g));

end