function [biglik,mu,sig,zm] = likEM_loop(inp,cons,expend,y,price,ot_ind,w_b,c_b,scl,v,type,zm)
%Calculates likelihood given type and weight of conspicuous consumption.  The coefficient on food at home is set to one, and all approximations are around one.

alp = inp(1);
mu = inp(2:29);
sig = inp(30:57);

%alp = exp(inp(1))/(1+exp(inp(1)));

%get the equilibrium c's for a wealth grid 
b = logical(1-eye(29));
g = @(x,K,k,m,bet) (1-alp)*sum(bet(b(:,k)))/(bet(k) + alp*(sum(bet(b(:,k)))))*price(m,k)*x + K*x.^(-bet(k)/(alp*(sum(bet(b(:,k)))))); 
g_p = @(x,K,k,m,bet) 1/alp*((1-alp)*price(m,k) - bet(k)/(sum(bet(b(:,k))))*g(x,K,k,m,bet)./x);

%Create constant and wealth/consumption grids (for ease of finding c given
%W)

K = cell(29,18);
for m = 1:18
    for k=1:29 
        K{k,m} = @(bet) (w_b(m)-(bet(k)+sum(bet(b(:,k))))/(bet(k) + alp*sum(bet(b(:,k))))*price(m,k)*bet(k)/sum(bet)*w_b(m)/price(m,k))*(bet(k)/sum(bet)*w_b(m)/price(m,k))^(bet(k)/(alp*sum(bet(b(:,k)))));
    end
end

%% calculate likelihood
%ot is observation type, should have one for each household
biglik = 0;
weirdo = 0;
ob_type = zeros(size(expend));
for h = 1:size(expend,1)
    lik = zeros(29,1);
    gamcd = cons(h,:)/cons(h,1);
    flag = zeros(29,1);
    for k = 1:29
        %display(h);
        ot = k;
        if ot ~= 1
            gam = gamcd;
            gam(ot) = sum(gam(b(:,ot)))*cons(h,ot)/price(y(h),ot)*((1-alp)*price(y(h),ot)/(expend(h)-cons(h,ot))-alp*g_p(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam)/g(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam));         
            if gam(ot)<0 || isnan(gam(ot)) == 1
                %weirdo = weirdo + 1;
                gam(ot) = gamcd(ot);
                flag(ot) = 1;
            end
                gl = gam(2:end);
                lik(ot) = log(v{type(h)}(ot))+sum((gl>0).*log(1-zm)) + sum((gl==0).*log(zm))+sum(log(max(lognpdf(gl(gl>0),mu(gl>0),sig(gl>0)),1e-12)));          
                if flag(ot) == 1
                    lik(k) = lik(k) - abs(lik(k));
                end
        else
%         gam = gamcd;
%         hat = (cons(h,ot)/price(y(h),ot)*((1-alp)*price(y(h),ot)/(expend(h)-cons(h,ot))-alp*g_p(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam)/g(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam)))^(-1);        
%         if hat < 0 || isnan(hat) == 1
%             %weirdo = weirdo + 1;
%             lik(k) = -100500; %punishment, large enough to hurt but not kill
%         else
%             gam(2:end) = cons(h,2:end)*hat/(expend(h)-cons(h,ot));
%             gl = gam(2:end);
%             if isempty(gl>0) == 1;
%                 lik(k) = log(v{type(h)}(k))+ sum((gl==0).*log(zm));
%             else
%                 lik(k) = log(v{type(h)}(k))+sum((gl>0).*log(1-zm)) + sum((gl==0).*log(zm))+sum(log(lognpdf(gl(gl>0),mu(gl>0),sig(gl>0))));
%             end
%         end
        end
    end
lik(1) = -inf; %suppose you can't be obs type 1
[addme,ob_type(h)] = max(lik);    
if min(flag(2:end)) == 1
    weirdo = weirdo + 1;
end
%if addme == -100000
%    weirdo = weirdo + 1;
%    ob_type(h) = -1;
%end
biglik = biglik +addme;  
end

% ngam = zeros(29,size(expend,1)); %gams,household
% for h = 1:size(expend,1)
%     %display(h);
%     ot = ot_ind(h);
%     if ot ~= 1
%         gam = cons(h,:)/cons(h,1);
%          gam(ot) = sum(gam(b(:,ot)))*cons(h,ot)/price(y(h),ot)*((1-alp)*price(y(h),ot)/(expend(h)-cons(h,ot))-alp*g_p(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam)/g(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam));         
%          if gam(ot)<0 || isnan(gam(ot)) == 1
%              flag(h) = 1;
%              weirdo = weirdo + 1;
%          end
%          ngam(:,h) = gam;
%     else
%         gam = cons(h,:)/cons(h,1);
%         hat = (cons(h,ot)/price(y(h),ot)*((1-alp)*price(y(h),ot)/(expend(h)-cons(h,ot))-alp*g_p(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam)/g(cons(h,ot)/price(y(h),ot),K{ot,y(h)}(gam),ot,y(h),gam)))^(-1);        
%         if hat < 0 || isnan(hat) == 1
%             flag(h) = 1;
%             weirdo = weirdo + 1;
%         end
%         gam(2:end) = cons(h,2:end)*hat/(expend(h)-cons(h,ot));
%         ngam(:,h) = gam;
%     end
% end
% 
% display(['A total of ', num2str(weirdo), ' weirdos.']);
% 
% mu = zeros(28,1);
% sig = zeros(28,1);
% zm = zeros(28,1);
% 
% lik = 0;
% for m = 2:29
%     k = m-1;
%     use = (flag==0)';
%     nmp = ngam(m,:)>0; %not mass points
%     mp = 1-nmp; %mass points
%     mu(k) = mean(log(ngam(m,nmp & use)));
%     sig(k) = sqrt(mean((log(ngam(m,nmp & use))-mu(k)).^2));
%     zm(k) = sum(mp & use)/sum(use);
%     lik = lik + sum(log(lognpdf(ngam(m,nmp & use),mu(k),sig(k))))+log(zm(k))*sum(mp & use)+log(1-zm(k))*sum(nmp & use);
% end
% 
% %Punish
% lik = lik - sum(flag)*100000;

biglik = -biglik/scl; %we want to maximize the likelihood, not minimize it!
ob_types = histc(ob_type(ob_type>0),1:29)/sum(ob_type>0);

display(mu);
display(sig);
display(ob_types);
display(biglik);
display(alp);
display(weirdo)
end

