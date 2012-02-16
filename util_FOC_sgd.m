function [u,foc] = util_FOC_sgd(p,cf,w,price,scl)
%function [u,foc,hessian] = util_FOC(p,input,g,alp,v,w,expend)
%calculates utility FOC for parameter vector p, approx coefficients cf and
%weath level w

cf = reshape(cf,size(cf,1)/28,28)';
for k = 1:size(cf,2)
    cf(:,k) = cf(:,k)/scl^(k-1);
end

u = 0;
grad = 0;

%wealth poly coefs
for j = 1:size(w,2)
    w_poly = zeros(size(cf,2),1);
    for k = 1:size(cf,2)
        w_poly(k) = w(j)^(k-1);
    end

    %get budget shares
    sh = zeros(29,1);
    for k = 1:28
        sh(k+1) = exp(cf(k,:)*w_poly);
    end
    sh = sh/(1+sum(sh));
    sh(1) = 1 - sum(sh);

    %get actual consumption
    cns = (w(j) * sh)./price';

    u = u + sum(bsxfun(@times,p(:,1),log(bsxfun(@plus,p(:,2),cns))),1);
    
    dudc = p(:,1)./bsxfun(@plus,p(:,2),cns);
    dcds = w(j)./price';
    sh_stack = repmat(sh',size(cf,2),1);
    dsdcf = -bsxfun(@times,bsxfun(@times,repmat(sh',size(cf(:),1),1),repmat(w_poly,28,1)),sh_stack(size(cf,2)+1:end)'); %size(#coefs,#goods)
    for k = 1:size(cf,2)
        dsdcf(size(dsdcf,1)+k:size(dsdcf,1)+size(cf,2):end) = dsdcf(size(dsdcf,1)+k:size(dsdcf,1)+size(cf,2):end)' + sh(2:end)*w_poly(k);
    end
    dudcf = sum(bsxfun(@times,dudc'.*dcds',dsdcf),2);
    
    grad = grad + dudcf(:);
    
end

%fundamental utility grad
%fu = (bsxfun(@times,p(:,1)',1./bsxfun(@plus,p(:,2)',c')))';

%equilibrium expenditure-utility vector
%guv = sum(bsxfun(@times,p(:,1)',log(bsxfun(@plus,p(:,2)',g))),2);

%likelihood vector for each wealth level
%dm = ones(size(g,1),size(g,2));
%dnsobj = bsxfun(@minus,g,c');
%jac_dm = ones(size(g,2),size(g,2),size(g,1));

%w_product = prod(w_prob,2);
%for m = 1:29
%   dm(:,m) = 1/(2*pi*v(m)^2)*(dnsobj(:,m)/(v(m)^2)).*exp(-dnsobj(:,m).^2/(2*v(m)^2));
%   dm(:,m) = dm(:,m).*w_product./w_prob(:,m); 
%end
%dens = bsxfun(@times,dm,w)/sum(w_product.*w);

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

%ex = sum(bsxfun(@times,guv,dens))'; %social expectation of utility level

%foc = fu;

% hessian = [[(1-alp)*bsxfun(@times,eye(29),-fu./bsxfun(@plus,p(:,2),c))+alp*sum(jac_dens,3);-ones(1,29)],[-ones(29,1);0]];
u = -u;
foc = -grad;
for k = 1:size(cf,2)
   foc(k:size(cf,2):end) = foc(k:size(cf,2):end)/scl^(k-1); 
end

end