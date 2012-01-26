function [u,foc] = util_FOC(p,c,g,alp,v,w)
%function [u,foc,hessian] = util_FOC(p,input,g,alp,v,w,expend)
%calculates utility FOC for parameter vector p, consumption vector c,
%equilibrium consumption matrix g:exp->cons, consp cons weight alp, and 
%shock std dev vector v, and weath density w

%c = exp(input(1:29));
c = c';
%sp = input(30);

[u,guv,w_prob,dnsobj] = util(p,c,g,alp,v,w);

%fundamental utility grad
fu = (bsxfun(@times,p(:,1)',1./bsxfun(@plus,p(:,2)',c')))';

%equilibrium expenditure-utility vector
%guv = sum(bsxfun(@times,p(:,1)',log(bsxfun(@plus,p(:,2)',g))),2);

%likelihood vector for each wealth level
dm = ones(size(g,1),size(g,2));
%dnsobj = bsxfun(@minus,g,c');
%jac_dm = ones(size(g,2),size(g,2),size(g,1));

w_product = prod(w_prob,2);
for m = 1:29
   dm(:,m) = -1/(2*pi*v(m)^2)*(dnsobj(:,m)/(v(m)^2)).*exp(-dnsobj(:,m).^2/(2*v(m)^2));
   dm(:,m) = dm(:,m).*w_product./w_prob(:,m); 
end
dens = bsxfun(@times,dm,w);

% for m = 1:29
%     for l = 1:m
%         k = 1;
%         while k<30 && k ~= m && k ~= l
%             temp = normpdf(dnsobj(:,k),0,v(k)).*temp;
%             k = k+1;
%         end
%         if m == l
%             jac_dm(m,l,:) =temp.*((-1/(2*pi*v(m)^2)*(1/(v(m)^2))*exp(-dnsobj(:,m).^2/(2*v(m)^2))+1/(2*pi*v(m)^2)*(dnsobj(:,m).^2/(v(m)^4)).*exp(-dnsobj(:,m).^2/(2*v(m)^2))));
%         else
%             temp = ones(size(g,1),1);
%             jac_dm(m,l,:) = prod([-1/(2*pi*v(m)^2)*(dnsobj(:,m)/(v(m)^2)).*exp(-dnsobj(:,m).^2/(2*v(m)^2)),...
%                 -1/(2*pi*v(l)^2)*(dnsobj(:,l)/(v(l)^2)).*exp(-dnsobj(:,l).^2/(2*v(l)^2)),temp],2);
%             jac_dm(l,m,:) = jac_dm(m,l,:);
%         end
%     end
% end
% jac_dens = zeros(size(jac_dm));
% for k = 1:29
%     temp = squeeze(jac_dm(k,:,:));
% jac_dens(k,:,:) = bsxfun(@times,temp,w'.*guv');
% end

ex = sum(bsxfun(@times,guv,dens))'; %social expectation of utility level

foc = (1-alp)*fu+alp*ex;

% hessian = [[(1-alp)*bsxfun(@times,eye(29),-fu./bsxfun(@plus,p(:,2),c))+alp*sum(jac_dens,3);-ones(1,29)],[-ones(29,1);0]];
u = -u;
foc = -foc;

end