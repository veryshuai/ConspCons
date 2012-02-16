function r = mpec_err(X,price,egr,v,w,ca,egr_cuts)

warning off all;

tic

%diary off
%diary 2-4-12.txt

param1 = X(1:29)'; %cobb douglas parameters
param2 = X(30:58)'; %stone geary parameters
alp = X(59); %weight of conspicuous consumption
vm = X(60); %v multiplier

%alp = 0;

v = v*vm; %scale v by the multiplier

%sobol points
s = sobolset(1);
s = net(s,10);

options=optimset('Display','iter','jacobian','off','TolFun',1e-2,'TolX',1e-2,'TolCon',1e-6,'DerivativeCheck','off',...
    'GradObj','off');%,'Algorithm','active-set');

init = ones(28*3,1);
eq = zeros(28*3,18,4);
for k = 1:18
    for m = 1:4
        [val,eq(:,k,m)] = fminunc(@(x) inner([param1,param2],x,w(:,k,m),price(k,:),alp,v(m,:),egr(k,:),s),init,options);
        init = eq(:,k,m);
    end
end

% display(param1);
% display(param2);
% display(alp);
% display(vm);



toc

end