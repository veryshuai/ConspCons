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
w_prob = ones(size(guv,1),29)*eps;
calc = ones(size(guv)); %indexes to calculate

%COMMENTED OUT: ATTEMPT TO INCREASE SPEED, BUT ACTUALLY FASTER THE ORIGINAL
%WAY
%ind = (1:size(guv,1))';
%false_ind = true(size(guv,1));
%tic
for k=1:29
    %calc(dnsobj(calc==1,k)>2*v(k) | dnsobj(calc==1,k)<-2*v(k))=0;
    w_prob(calc==1,k) = normpdf(dnsobj(calc==1,k),0,v(k)); 
    % ind = find(w_prob(:,k)>1e-5);
    % false_ind(ind) = false;
    % w_prob(false_ind,k) = 0;
   dm = w_prob(:,k).*dm;
end
den_u = sum(dm.*w);
dens = dm.*w/den_u; %make it a density
dens = max(eps,dens); %make sure there are no zeros

ex = sum(guv.*dens); %social expectation of utility level

u = (1-alp)*fu+alp*ex;
%display(u);
%toc
end