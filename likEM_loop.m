function lik = likEM_loop(alp,cons,exp,y,price,ot_ind,w_b,c_b,scl)
%Calculates likelihood given type and weight of conspicuous consumption.  The coefficient on food at home is set to one, and all approximations are around one.

%get the equilibrium c's for a wealth grid 
b = logical(1-eye(29));
g = @(x,K,k,m,bet) (1-alp)*(sum(bet(b(:,k)))-bet(k))/(bet(k) + alp*(sum(bet(b(:,k)))-bet(k)))*x + K*x.^(bet(k)/(alp*(sum(bet(b(:,k)))-bet(k)))); 
g_w =  @(x,K,k,m,bet) g(x,K,k,m,bet) + price(m,k)*x;
g_p = @(x,K,k,m,bet) 1/alp*((1-alp)*price(m,k) - bet(k)/(sum(bet(b(:,k))))*g(x,K,k,m,bet)./x);
g_pp = @(x,K,k,m,bet) -bet(k)/(alp*(sum(bet(b(:,k)))))*(g_p(x,K,k,m,bet)./x-g(x,K,k,m,bet)./(x.^2));

%Create constant and wealth/consumption grids (for ease of finding c given
%W)

K = cell(29,18);
cg = cell(29,18);
wg = cell(29,18);
for m = 1:18
    for k=1:29 
        %K{k,m} = @(bet) (w_b(m)-((sum(bet(b(:,k)))-bet(k))*(1-alp)/(bet(k)+alp*(sum(bet(b(:,k)))-bet(k)))+price(m,k))*bet(k)/price(m,k)*w_b(m))/(bet(k)/price(m,k)*w_b(m))^(bet(k)/((sum(bet(b(:,k)))-bet(k))*alp));   
        K{k,m} = @(bet) (bet(k)/(price(m,k)*sum(bet(b(:,k))))*w_b(m))^(1-bet(k)/(alp*sum(bet(b(:,k)))))*(1-sum(bet(b(:,k)))*(1-alp)/(bet(k)+alp*sum(bet(b(:,k)))));
        cg{k,m} = linspace(1e-12,c_b(m),1000);
        wg{k,m} = @(bet) g_w(cg{k,m},K{k,m}(bet),k,m,bet);
    end
end

%create dgdc as a function of other gammas, for not food at home
dgdc = cell(29,18);
for m = 1:18
    for k=1:29     
        dgdc{k,m} = @(x,gam,W,k,m) gam(k)./x + alp * x * sum(gam(b(:,k))) * g_p(x,K{k,m}(gam),k,m,gam).^2/g(x,K{k,m}(gam),k,m,gam).^2 ...
            -alp * x * sum(gam(b(:,k))) * g_pp(x,K{k,m}(gam),k,m,gam)/g(x,K{k,m}(gam),k,m,gam) ...
            +price(m,k).^2 * x ./((W - price(m,k) * x).^2 );
    end
end

%create dhdc for not food at home
dhdc = cell(29,18);
b = logical(1-eye(29));
for m = 1:18
    for k=1:29     
        dhdc{k,m} = @(x,gam,W,k,m) (gam(k)./x + price(m,k).^2 ./(W - price(m,k) * x)) * g(x,K{k,m}(gam),k,m,gam)./(alp * g_p(x,K{k,m}(gam),k,m,gam)) + sum(gam(b(:,k))) * (g_p(x,K{k,m}(gam),k,m,gam)/g(x,K{k,m}(gam),k,m,gam));% - g_pp(x,K{k,m}(gam),k,m,gam)/g_p(x,K{k,m}(gam),k,m,gam));
    end
end

%% calculate likelihood
%ot is observation type, should have one for each household
%tic
nan_count = 0;
inf_count = 0;
sv = zeros(29,1); %best sample variances
ngam = zeros(29,size(exp,1)); %gams,household
for h = 1:size(exp,1)
    %display(h);
    ot = ot_ind(h);
    if ot ~= 1
        gam = cons(h,:)/cons(h,1);
        gam(ot) = 1;
        tempw = wg{ot,y(h)}(gam);
    	ind = find(exp(h)>tempw,1,'last');
    	ogam = 1 + dgdc{ot,y(h)}(cg{ot,y(h)}(ind),gam,tempw(ind),ot,y(h))*(cons(h,ot)/price(y(h),ot)-cg{ot,y(h)}(ind));
        gam(ot) = max(ogam,0);
        ngam(:,h) = gam;
    else
        gam = ones(29,1);
        tempw = wg{ot,y(h)}(gam);
        ind = find(exp(h)>tempw,1,'last');
        hgam = 28 + dhdc{ot,y(h)}(cg{ot,y(h)}(ind),gam,tempw(ind),ot,y(h))*(cons(h,ot)/price(y(h),ot)-cg{ot,y(h)}(ind));
        hat = max(hgam,0);
        gam(2:end) = cons(h,2:end)*hat/(exp(h)-cons(h,ot));
        ngam(:,h) = gam;
    end
end

%to deal with very unlikely guys...
ngam = min(ngam,100);

lam = zeros(29,1);
lik = 0;
for k = 1:29
    lam(k) = 1/mean(ngam(k,:));
    lik    = lik + sum(log(exppdf(ngam(k,:),lam(k))));
end

%%likmin = 0;
%nan_count = 0;
%inf_count = 0;
    
%    temp = sum(log(normpdf(ngam(:,h),zeros(29,1),sv)));
%    if isnan(temp)==0 && abs(temp)<inf
%        lik = lik + temp;
%        %likmin = min(temp,likmin);
%    elseif isnan(temp) == 0
%        nan_count = nan_count + 1;
%        lik = lik + -1e4; %don't want to give an incentive to create infs and nans
%    else
%        inf_count = inf_count + 1;
%        lik = lik + -1e4; %don't want to give an incentive to create infs and nans
%    end
%end
%
lik = -lik/scl; %we want to maximize the likelihood, not minimize it!
%display(nan_count);
%display(inf_count);
%toc
end
