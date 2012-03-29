function not = not_ind(X,cons,exp,y,type,price,w_b,c_b,v)
%Calculates likelihood given type

param1 = X(1:29); %cobb douglas parameters
param2 = X(30:58); %stone geary parameters
alp = X(59); %weight of conspicuous consumption
%sig = X(60:88); %std devs of preference shocks

%% prepare for solving for the invisible goods
bet = param1/sum(param1);
slv = cell(18,2);
for k = 1:18
 slv{k,1} = -bet*price(k,:);
 slv{k,1}(1:30:end) = (1-bet').*price(k,:);
 slv{k,2} = @(x,y,z) price(k,:).*x - price(k,:).*param2' - bet'*y + bet'*(price(k,:)*param2)+bet'*(price(k,z)*x(z));
end

%% prepare for solving visible goods

%get the equilibrium c's for a wealth grid 

g = @(x,K,k,m) price(m,:)*param2 - param2(k)*price(m,k) + (1-alp)*(1-bet(k))/(bet(k) + alp*(1-bet(k)))*(x-param2(k)) + K*(x-param2(k)).^(bet(k)/(alp*(1-bet(k)))); 
g_w =  @(x,K,k,m) g(x,K,k,m) + price(m,k)*x;
g_p = @(x,K,k,m) 1/alp*((1-alp)*price(m,k) - bet(k)/(1-bet(k))*(g(x,K,k,m)-price(m,:)*param2 + price(m,k)*param2(k))./(x-param2(k)));
g_pp = @(x,K,k,m) bet(k)/(alp*(1-bet(k))*(x-param2(k)).^2)*(g(x,K,k,m)-price(m,:)*param2 + price(m,k)*param2(k)-g_p(x,K,k,m)*(x-param2(k)));

%Create constant and wealth/consumption grids (for ease of finding c given
%W)

K = zeros(29,18);
cg = cell(29,18);
wg = cell(29,18);
for m = 1:18
    for k=1:29     
        K(k,m) = (w_b(m)-price(m,:)*param2-(1-bet(k))*(1-alp)/(bet(k)+alp*(1-bet(k)))*bet(k)/price(m,k)*(w_b(m)-price(m,:)*param2))/(bet(k)/price(m,k)*(w_b(m)-price(m,:)*param2))^(bet(k)/((1-bet(k))*alp));   
        cg{k,m} = linspace(param2(k),c_b(m),1000);
        wg{k,m} = g_w(cg{k,m},K(k,m),k,m);
    end
end

%create dpdc as a function of other gammas
dpdc = cell(29,18);
b = logical(1-eye(29));
for m = 1:18
    for k=1:29     
        dpdc{k,m} = @(x,gam,W,k,m) 1 + alp * (1-bet(k))/bet(k) *(x-param2(k)).^1 * (-g_p(x,K(k,m),k,m).^2/(g(x,K(k,m),k,m)...
            -price(m,b(:,k))*param2(b(:,k))).^2+g_pp(x,K(k,m),k,m)/(g(x,K(k,m),k,m)-price(m,:)*param2-price(m,k)-param2(k)) ...
            +price(m,k)/((W - price(m,b(:,k))*param2(b(:,k)) - price(m,b(:,k))*gam - price(m,k) * x).^2 * alp * (1-bet(k))));
    end
end

%% calculate likelihood
%ot is observation type, should have one for each household
not = zeros(size(exp,1),1);
sv = zeros(29); %rows are best sample variances, columns are observation types 
nan_count = zeros(29,1);
inf_count = zeros(29,1);
ngam = zeros(29,size(exp,1),29); %gams,household,ob_type
for h = 1:size(exp,1)
    %display(h)
    for k = 1:29
        ot = k;
        temp = slv{y(h),2}(cons(h,:),exp(h),ot);
        gam = slv{y(h),1}(b(:,ot),b(:,ot))\temp(b(:,ot))'; 
        ind = find(exp(h)>wg{ot,y(h)},1,'last');
        ogam = dpdc{ot,y(h)}(cg{ot,y(h)}(ind),gam,wg{ot,y(h)}(ind),ot,y(h))*(cons(h,ot)-cg{ot,y(h)}(ind));
        ins = 1:28;
        ins_ind = ins<ot;
        ngam(:,h,k) = [gam(ins_ind); ogam; gam(~ins_ind)];
        temp = ngam(:,h,k).^2;
        if max(isnan(temp)) == 0 && max(abs(temp)) ~= inf
           sv(:,k) = sv(:,k) + temp;
        elseif isnan(temp) == 1
           nan_count(k) = nan_count(k) + 1; 
        else
           inf_count(k) = inf_count(k) + 1;
        end
    end
end

for k = 1:29
    sv(:,k) = sqrt(sv(:,k)/(size(exp,1)-inf_count(k)-nan_count(k)));
end

parfor h = 1:size(exp,1)
    lik = zeros(29,1);
    for k = 1:29
        temp = sum(log(normpdf(ngam(:,h,k),zeros(29,1),sv(:,k))));
        if isnan(temp)==0 && abs(temp)<inf
            lik(k) = temp*v{type(h)}(k);
        elseif isnan(temp) == 1
            lik(k) = -inf;
%            nan_count = nan_count + 1;
        else
            lik(k) = -inf;
%            inf_count = inf_count + 1;
        end
    end
    [~,ind] = max(lik);
    not(h) = ind;
end
%display(nan_count);
%display(inf_count);

end