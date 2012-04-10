%This script runs welfare comparisons for conspicuous consumption vs no
%conspicuous consumption

clear;
load 3-23-12.mat;

param1 = init(1:29); %cobb douglas parameters
param2 = init(30:58); %stone geary parameters
alp = init(59); %weight of conspicuous consumption
%sig = X(60:88); %std devs of reference shocks

param1 = param1/sum(param1);


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


%wellev = figure(1);

%plot(regwel(1:50)-ccwel(1:50),'r');

cvec = {'y' 'm' 'c' 'r' 'g' 'b'};

for k = 1:21
m = 16;
%Welfare with consp cons
bh = sum(param1) - param1(k);
gh = sum(price(m,:)*param2-price(m,k)*param2(k));
z = log(prod((param1./price(m,:)').^param1)/((param1(k)/price(m,k))^param1(k))*(bh)^(-bh));

ccwel = bh*log(wg{k,m}-price(m,k)*cg{k,m}-gh)+param1(k)*log(cg{k,m}-param2(k)) + z;
regwel = 0;

for n = 1:29
   regwel = regwel + param1(n)*log(param1(n)/(sum(param1)*price(m,n))*(wg{k,m}-price(m,:)*param2)); 
end

wellos = figure(2);

plot((regwel(1:50)-ccwel(1:50))./regwel(1:50),cvec{mod(k,6)+1})

hold on

end

hold off
set(gca,'XTickLabel','');
xlabel('Expenditure');
ylabel('Percentage Utility Loss');
title('Utility Loss due to Conspicuous Consumption (by observtion type)')
print(wellos,'-dpdf','figures/wellos');
