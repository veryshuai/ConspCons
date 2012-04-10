function [sig_util,no_sig_util] = util_calc(wealth,alp,price,avp,k)
%function [mutil,util] = util_calc(tax,price,gam,w,ot,alp)
%this function calculates utility for a vector of wealth levels 
display(k);
avp = [1;avp];

b = logical(1-eye(29));
g = @(x,K,k,bet) (1-alp)*sum(bet(b(:,k)))/(bet(k) + alp*(sum(bet(b(:,k)))))*price(k)*x + K*x.^(-bet(k)/(alp*(sum(bet(b(:,k)))))); 
g_p = @(x,K,k,bet) 1/alp*((1-alp)*price(k) - bet(k)/(sum(bet(b(:,k))))*g(x,K,k,bet)./x);
K = @(k,bet) (1000-(bet(k)+sum(bet(b(:,k))))/(bet(k) + alp*sum(bet(b(:,k))))*price(k)*bet(k)/sum(bet)*1000/price(k))*(bet(k)/sum(bet)*1000/price(k))^(bet(k)/(alp*sum(bet(b(:,k)))));

foc = @(x,k,gam,W) (-(1-alp)*sum(gam(b(:,k)))/(W-price(k)*x) +alp*g_p(x,K(k,gam),k,gam)/g(x,K(k,gam),k,gam) + gam(k)/x)^2;

options = optimset('Display','off','Algorithm','active-set');
c_sig = zeros(size(wealth,1),29);
c_no_sig = zeros(size(wealth,1),29);
sig_util = zeros(size(wealth,1),1);
no_sig_util = zeros(size(wealth,1),1);
%init = max(avp(k)/sum(avp)/price(k),.01);
init = 1;
for w = 1:size(wealth,1)
    %display(w);
    c_sig(w,k) = fmincon(@(x) foc(x,k,avp,wealth(w)),init,[],[],[],[],0,wealth(w)/price(k),[],options); 
    %display(c_sig(w,k));
    %init = c_sig(w,k);
    c_no_sig(w,:) = avp'/sum(avp)./price * wealth(w); 
    index = 1:29;
    below = index<k;
    above = index>k;            
    not = logical(below + above);
    if wealth(w) >= c_sig(w,k)*price(k) && isnan(c_sig(w,k)) == 0
        c_sig(w,not) = (wealth(w)-c_sig(w,k)*price(k)) * avp(not)'/sum(avp(not))./price(not);
        sig_util(w) = avp(c_sig(w,:)>0)'*log(c_sig(w,c_sig(w,:)>0))';
    else 
        c_sig(w,:) = zeros(1,29);
        sig_util(w) = NaN;
        %flag(k) = 1; 
    end
    no_sig_util(w) = avp(c_no_sig(w,:)>0)'*log(c_no_sig(w,c_no_sig(w,:)>0))';
end

end

