function [mutil,util,adjc] = ag_util(tax,price,gam,w,ot,alp)
%this function calcuates aggregate utility for a given tax 

%oprice = price;
price = price .* tax';

b = logical(1-eye(29));
g = @(x,K,k,bet) (1-alp)*sum(bet(b(:,k)))/(bet(k) + alp*(sum(bet(b(:,k)))))*price(k)*x + K*x.^(-bet(k)/(alp*(sum(bet(b(:,k)))))); 
g_p = @(x,K,k,bet) 1/alp*((1-alp)*price(k) - bet(k)/(sum(bet(b(:,k))))*g(x,K,k,bet)./x);
K = @(k,bet) (1000-(bet(k)+sum(bet(b(:,k))))/(bet(k) + alp*sum(bet(b(:,k))))*price(k)*bet(k)/sum(bet)*1000/price(k))*(bet(k)/sum(bet)*1000/price(k))^(bet(k)/(alp*sum(bet(b(:,k)))));

foc = @(x,k,gam,W) (-(1-alp)*sum(gam(b(:,ot(k))))/(W-price(ot(k))*x) +alp*g_p(x,K(k,gam),k,gam)/g(x,K(k,gam),k,gam) + gam(ot(k))/x)^2;

c = zeros(size(ot,1),29);
util = zeros(size(ot,1),1);
options = optimset('Display','off','Algorithm','active-set');
flag = zeros(size(ot,1),1);
for k = 1:size(ot,1)
    %first get consumption of the observed good
    %init = max(gam(k,ot(k))/sum(gam(k,:))/price(ot(k)),1);
    init = 1;
    c(k,ot(k)) = fmincon(@(x) foc(x,ot(k),gam(k,:),w(k)),init,[],[],[],[],0,w(k)/price(ot(k)),[],options); 
    index = 1:29;
    below = index<ot(k);
    above = index>ot(k);            
    not = logical(below + above);
    if w(k) >= c(k,ot(k))*price(ot(k)) && isnan(c(k,ot(k))) == 0
        c(k,not) = (w(k)-c(k,ot(k))*price(ot(k))) * gam(k,not)/sum(gam(k,not));
        %util(k) = gam(c(k,:)>0)*log(c(k,c(k,:)>0))';
    else 
        c(k,:) = zeros(1,29);
        flag(k) = 1; 
    end
end

% Calculate Tax revenues
agc = sum(c,1);
rev = sum((tax'-1) .* agc);
inc = rev/sum(agc);
adjc = c /(1-inc); %geometric sum...cool!

%calculate mean utility (so I don't have to worry about flags!)
for k = 1:size(ot,1)
    if flag(k) == 0
        util(k) = gam(adjc(k,:)>0)*log(adjc(k,adjc(k,:)>0)./price(adjc(k,:)>0))';
    end
end
%display(sum(flag));
mutil = -sum(util)/(size(ot,1)-sum(flag));
end

