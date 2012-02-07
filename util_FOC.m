function [u,foc] = util_FOC(p,cf,w,price,scl,gcf,alp,s,v,egr)
%function [u,foc,hessian] = util_FOC(p,input,g,alp,v,w,expend)
%calculates utility FOC for parameter vector p, approx coefficients cf and
%weath level w


%% get sobol stuff  (this should ultimately be outside of this loop)
gcf = reshape(gcf,size(gcf,1)/28,28)';
for k = 1:size(gcf,2)
    gcf(:,k) = gcf(:,k)/scl^(k-1);
end

g  = zeros(size(s,1),29);
su = zeros(size(s,1),1);
ws = zeros(size(s,1),1);

for j = 1:size(s,1)
    
    w_poly = zeros(size(gcf,2),1);
    for k = 1:size(gcf,2)
        w_poly(k) = s(j)^(k-1);
    end
    
    %density
    ws(j) = w_poly'*w(1:end-1);%/w(end);
    
    %get budget shares
    sh = zeros(29,1);
    for k = 1:28
        sh(k+1) = exp(gcf(k,:)*w_poly);
    end
    sh = sh/(1+sum(sh));
    sh(1) = 1 - sum(sh);

    %get actual consumption
    g(j,:) = (s(j) * sh)./price';
    
    %get sobol utilities
    su(j) = sum(bsxfun(@times,p(:,1),log(bsxfun(@plus,p(:,2),g(j,:)'))),1);
end

%% get utility level for each wealth type 

cf = reshape(cf,size(cf,1)/28,28)';
for k = 1:size(cf,2)
    cf(:,k) = cf(:,k)/scl^(k-1);
end

tu = 0;
tg = 0;

%wealth poly coefs
for j = 1:size(egr,2)
    w_poly = zeros(size(cf,2),1);
    for k = 1:size(cf,2)
        w_poly(k) = egr(j)^(k-1);
    end

    %get budget shares
    sh = zeros(29,1);
    for k = 1:28
        sh(k+1) = exp(cf(k,:)*w_poly);
    end
    sh = sh/(1+sum(sh));
    sh(1) = 1 - sum(sh);

    %get actual consumption
    cns = (egr(j) * sh)./price';

    fu = sum(bsxfun(@times,p(:,1),log(bsxfun(@plus,p(:,2),cns))),1);
    
    diff = bsxfun(@minus,g,cns');
    wp_int_int = normpdf(diff,zeros(size(diff)),repmat(v,size(diff,1),1));
    wp_int = bsxfun(@times,wp_int_int,ws);
    wp = prod(wp_int,2)/sum(prod(wp_int,2));
    ex = sum(su.*wp);
    
    tu = tu + (1-alp)*fu + alp * ex;
    
    %get fundy grad
    dudc = p(:,1)./bsxfun(@plus,p(:,2),cns);
    dcds = egr(j)./price';
    sh_stack = repmat(sh',size(cf,2),1);
    dsdcf = -bsxfun(@times,bsxfun(@times,repmat(sh',size(cf(:),1),1),repmat(w_poly,28,1)),sh_stack(size(cf,2)+1:end)'); %size(#coefs,#goods)
    for k = 1:size(cf,2)
        dsdcf(size(dsdcf,1)+k:size(dsdcf,1)+size(cf,2):end) = dsdcf(size(dsdcf,1)+k:size(dsdcf,1)+size(cf,2):end)' + sh(2:end)*w_poly(k);
    end
    dcdcf = bsxfun(@times,dcds',dsdcf);
    dudcf = sum(bsxfun(@times,dudc',dcdcf),2);
    grad = dudcf;
    
    %get integral grad
    dddc = zeros(size(s,1),29,29);
    for k = 1:size(s,1)
        for m = 1:29
        dddc(:,m,m) = 1;
        end
    end
    
    v_rep = repmat(v,size(s,1),1);
    njacs = 1./(2*pi*v_rep.^2).*(diff./(v_rep.^2)).*exp(-diff.^2./(2*v_rep.^2));
    fot_int = zeros(size(s,1),29);
    fot_int(:,1) = prod(bsxfun(@times,wp_int_int(:,2:end),ws),2);
    for k = 2:28
       fot_int(:,k) = prod(bsxfun(@times,wp_int_int(:,[1:k-1,k+1:end]),ws),2);  
    end
    fot_int(:,29) = prod(bsxfun(@times,wp_int_int(:,1:end-1),ws),2);
    fot = njacs.*fot_int/sum(prod(wp_int,2));
    sot_int = repmat(-wp/sum(prod(wp_int,2)),[1,29]);
    sot = sot_int.*fot_int.*njacs;
    dwdd = fot + sot;
    
    dedw = su;
    
    dedc = zeros(29,1);
    for k = 1:29
       dedc(k) = dedw'*sum(dwdd.*dddc(:,:,k),2); 
    end
    
    dedcf = sum(bsxfun(@times,dedc',dcdcf),2);
    
    grad_ex = dedcf;
    
    tg = tg + (1-alp)*grad + alp*grad_ex;
         %(:,m) = 1/(2*pi*v(m)^2)*(diff(:,m)/(v(m)^2)).*exp(-dnsobj(:,m).^2/(2*v(m)^2));
         %dm(:,m) = dm(:,m).*w_product./w_prob(:,m); 
%     end
%     dwdd = bsxfun(@times,dm,w)/sum(w_product.*w);
%     
%     dedw = 0;
%     
%     
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
u = - tu;
foc = -tg;
for k = 1:size(cf,2)
   foc(k:size(cf,2):end) = foc(k:size(cf,2):end)/scl^(k-1); 
end

end