function [lik,mu,sig,zm] = likEM_loop(inp,cons,expend,y,price,ot_ind,w_b,c_b,scl,v)
%Calculates likelihood given type and weight of conspicuous consumption.  The coefficient on food at home is set to one, and all approximations are around one.

alp = inp;
%alp = exp(inp(1))/(1+exp(inp(1)));
%zm  = exp(inp(2:end))./(1+exp(inp(2:end)));

%get the equilibrium c's for a wealth grid 
b = logical(1-eye(29));
% g = @(x,K,k,m,bet) (1-alp)*sum(bet(b(:,k)))/(bet(k) + alp*(sum(bet(b(:,k)))))*x + K*x.^(bet(k)/(alp*(sum(bet(b(:,k)))))); 
% g_w =  @(x,K,k,m,bet) g(x,K,k,m,bet) + price(m,k)*x;
% g_p = @(x,K,k,m,bet) 1/alp*((1-alp)*price(m,k) - bet(k)/(sum(bet(b(:,k))))*g(x,K,k,m,bet)./x);
% g_pp = @(x,K,k,m,bet) -bet(k)/(alp*(sum(bet(b(:,k)))))*(g_p(x,K,k,m,bet)./x-g(x,K,k,m,bet)./(x.^2));
% 
% %Create constant and wealth/consumption grids (for ease of finding c given
% %W)
% 
% K = cell(29,18);
% cg = cell(29,18);
% wg = cell(29,18);
% for m = 1:18
%     for k=1:29 
%         %K{k,m} = @(bet) (w_b(m)-((sum(bet(b(:,k))))*(1-alp)/(bet(k)+alp*(sum(bet(b(:,k))))+price(m,k))*bet(k)/price(m,k)*w_b(m))/(bet(k)/price(m,k)*w_b(m))^(bet(k)/((sum(bet(b(:,k))))*alp));   
%         K{k,m} = @(bet) (bet(k)/(price(m,k)*sum(bet))*w_b(m))^(1-bet(k)/(alp*sum(bet(b(:,k)))))*(1-sum(bet(b(:,k)))*(1-alp)/(bet(k)+alp*sum(bet(b(:,k)))));
%         cg{k,m} = linspace(1e-12,c_b(m)/price(m,k),1000);
%         wg{k,m} = @(bet) g_w(cg{k,m}/price(m,k),K{k,m}(bet),k,m,bet);
%     end
% end
% 
% %create dgdc as a function of other gammas, for not food at home
% dgdc = cell(29,18);
% for m = 1:18
%     for k=1:29     
%         %dgdc{k,m} = @(x,gam,W,k,m) gam(k)./x + alp * x * sum(gam(b(:,k))) * g_p(x,K{k,m}(gam),k,m,gam).^2/g(x,K{k,m}(gam),k,m,gam).^2 ...
%         %    -alp * x * sum(gam(b(:,k))) * g_pp(x,K{k,m}(gam),k,m,gam)/g(x,K{k,m}(gam),k,m,gam) ...
%         %    +price(m,k).^2 * x ./((W - price(m,k) * x).^2 );
%         dgdc{k,m} = @(x,gam,W,k,m) (1/x)/(gam(k)/x^2 - alp * sum(gam(b(:,k)))*g_pp(x,K{k,m}(gam),k,m,gam)/g(x,K{k,m}(gam),k,m,gam)...
%             + alp * sum(gam(b(:,k)))*g_p(x,K{k,m}(gam),k,m,gam)^2/g(x,K{k,m}(gam),k,m,gam)^2 + price(m,k)^2/(W-price(m,k)*x)^2);
%     end
% end
% 
% %create dhdc for not food at home
% dhdc = cell(29,18);
% b = logical(1-eye(29));
% for m = 1:18
%     for k=1:29     
%         dhdc{k,m} = @(x,gam,W,k,m) (gam(k)./x + price(m,k).^2 ./g(x,K{k,m}(gam),k,m,gam).^2) * g(x,K{k,m}(gam),k,m,gam)./(alp * g_p(x,K{k,m}(gam),k,m,gam)) + sum(gam(b(:,k))) * (g_p(x,K{k,m}(gam),k,m,gam)/g(x,K{k,m}(gam),k,m,gam) - g_pp(x,K{k,m}(gam),k,m,gam)/g_p(x,K{k,m}(gam),k,m,gam));
%     end
% end

%% calculate likelihood
%ot is observation type, should have one for each household
%tic
%nan_count = 0;
%inf_count = 0;
%sv = zeros(29,1); %best sample variances
ngam = zeros(29,size(expend,1)); %gams,household
for h = 1:size(expend,1)
    %display(h);
    ot = ot_ind(h);
    if ot ~= 1
        gam = cons(h,:)/cons(h,1);
        %gam(ot) = 1;
        %tempw = wg{ot,y(h)}(gam);
    	%ind = find(expend(h)>tempw,1,'last');
        %ogam = gam(ot) + dgdc{ot,y(h)}(cg{ot,y(h)}(ind)/price(y(h),ot),gam,tempw(ind),ot,y(h))*(cons(h,ot)-cg{ot,y(h)}(ind))/price(y(h),ot);
        %gam(ot) = max(ogam,0);
        gam(ot) = sum(gam(b(:,ot)))*((1-alp)/(expend(h)/(cons(h,ot)*price(y(h),ot))-1)-alp);
        if gam(ot)<0
            gam(ot) = 50; %punishment, give a large enough value to hurt, but not inf which kills likelihood
        end
        ngam(:,h) = gam;
    else
        gam = cons(h,:)/cons(h,1);
        %tempw = wg{ot,y(h)}(gam);
        %ind = find(expend(h)>tempw,1,'last');
        %hgam = sum(gam(2:end)) + dhdc{ot,y(h)}(cg{ot,y(h)}(ind)/price(y(h),ot),gam,tempw(ind),ot,y(h))*(cons(h,ot)-cg{ot,y(h)}(ind))/price(y(h),ot);
        %hat = max(hgam,0);
        hat = ((1-alp)/(expend(h)/(cons(h,ot)*price(y(h),ot))-1)-alp)^-1;
        if hat < 0
            hat = 500; %punishment, large enough to hurt but not kill
        end
        gam(2:end) = cons(h,2:end)*hat/(expend(h)-cons(h,ot));
        ngam(:,h) = gam;
    end
end

%to deal with very unlikely guys...
%ngam = min(ngam,100);

mu = zeros(28,1);
sig = zeros(28,1);
zm = zeros(28,1);

lik = 0;
for m = 2:29
    k = m-1;
    mu(k) = mean(log(ngam(m,ngam(m,:)>0)));
    sig(k) = sqrt(mean((log(ngam(m,ngam(m,:)>0))-mu(k)).^2));
    zm(k) = sum(ngam(m,:)==0)/size(ngam(m,:),2);
    lik = lik + sum(log(lognpdf(ngam(m,ngam(m,:)>0),mu(k),sig(k))))+log(zm(k))*sum(ngam(m,:)==0)+log(1-zm(k))*sum(ngam(m,:)>0);
    %lam(k) = 1/mean(ngam(k,:));
    %lik    = lik + sum(log(exppdf(ngam(k,:),lam(k))));
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
