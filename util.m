function [u,guv,w_prob,dnsobj] = util(p,c,g,alp,v,w)
%calculates utility value for parameter vector p, consumption vector c,
%equilibrium consumption matrix g:exp->cons, consp cons weight alp, and 
%shock std dev vector v, and weath density w

%fundamental utility
fu = sum(bsxfun(@times,p(:,1)',log(bsxfun(@plus,p(:,2)',c'))),2);

%equilibrium expenditure-utility vector
guv = sum(bsxfun(@times,p(:,1)',log(bsxfun(@plus,p(:,2)',g))),2);

%likelihood vector for each wealth level
dnsobj = bsxfun(@minus,g,c');
dm = ones(size(g,1),1);
w_prob = zeros(size(guv,1),29);
for k=1:29
   w_prob(:,k) = normpdf(dnsobj(:,k),0,v(k)); 
   dm = w_prob(:,k).*dm;
end
den_u = sum(dm.*w);
dens = dm.*w/den_u; %make it a density


ex = sum(guv.*dens); %social expectation of utility level

u = (1-alp)*fu+alp*ex;

end