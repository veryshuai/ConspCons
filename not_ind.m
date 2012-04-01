function not = not_ind(inp,cons,expend,y,type,price,w_b,c_b,v,mu,sig,zm)
%Calculates likelihood given type%Calculates likelihood given type and weight of conspicuous consumption.  The coefficient on food at home is set to one, and all approximations are around one.

alp = exp(inp)/(1+exp(inp));
%alp = inp;

%get the equilibrium c's for a wealth grid 
b = logical(1-eye(29));
g = @(x,K,k,m,bet) (1-alp)*(sum(bet(b(:,k))))/(bet(k) + alp*(sum(bet(b(:,k)))))*x + K*x.^(bet(k)/(alp*(sum(bet(b(:,k)))-bet(k)))); 
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
        K{k,m} = @(bet) (bet(k)/(price(m,k)*sum(bet))*w_b(m))^(1-bet(k)/(alp*sum(bet(b(:,k)))))*(1-sum(bet(b(:,k)))*(1-alp)/(bet(k)+alp*sum(bet(b(:,k)))));
        cg{k,m} = linspace(1e-12,c_b(m)/price(m,k),1000);
        wg{k,m} = @(bet) g_w(cg{k,m}/price(m,k),K{k,m}(bet),k,m,bet);
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
not = zeros(size(expend,1),1);
for h = 1:size(expend,1)
    lik = zeros(29,1);
    for k = 1:29
    %display(h);
    ot = k;
    if ot ~= 1
        gam = cons(h,:)/cons(h,1);
        %gam(ot) = 1;
        tempw = wg{ot,y(h)}(gam);
    	ind = find(expend(h)>tempw,1,'last');
    	ogam = gam(ot) + dgdc{ot,y(h)}(cg{ot,y(h)}(ind)/price(y(h),ot),gam,tempw(ind),ot,y(h))*(cons(h,ot)-cg{ot,y(h)}(ind))/price(y(h),ot);
        gam(ot) = max(ogam,0);
        gl = gam(2:end);
        lik(k) = log(v{type(h)}(k))+sum((gl>0).*log(1-zm')) + sum((gl==0).*log(zm'))+sum(log(lognpdf(gl(gl>0)',mu(gl>0),sig(gl>0))));
    else
        gam = cons(h,:)/cons(h,1);
        tempw = wg{ot,y(h)}(gam);
        ind = find(expend(h)>tempw,1,'last');
        hgam = sum(gam(2:end)) + dhdc{ot,y(h)}(cg{ot,y(h)}(ind)/price(y(h),ot),gam,tempw(ind),ot,y(h))*(cons(h,ot)-cg{ot,y(h)}(ind))/price(y(h),ot);
        hat = max(hgam,0);
        %if hat == 0;
        %    lik = -inf;
        %else
        gam(2:end) = cons(h,2:end)*hat/(expend(h)-cons(h,ot));
        gl = gam(2:end);
        if isempty(gl>0) == 1;
            lik(k) = log(v{type(h)}(k))+ sum((gl==0).*log(zm));
        else
            lik(k) = log(v{type(h)}(k))+sum((gl>0).*log(1-zm')) + sum((gl==0).*log(zm'))+sum(log(lognpdf(gl(gl>0)',mu(gl>0),sig(gl>0))));
        end
    end
    end
[~,not(h)] = max(lik); 
end


end