function [u,foc] = util_FOC(p,cf,w,price,scl,gcf,alp,s1,s2,v,egr,pscale)
%function [u,foc,hessian] = util_FOC(p,input,g,alp,v,w,expend)
%calculates utility FOC for parameter vector p, approx coefficients cf and
%weath level w

s = cell(2,1);
s{1} = s1;
s{2} = s2;
%% get numerical integration utilities (su, sobel utility)
gcf = reshape(gcf,size(gcf,1)/28,28)';
for k = 1:size(gcf,2)
    gcf(:,k) = gcf(:,k)/scl^(k-1);
end

g  = zeros(size(s{1},1),29);
su = zeros(size(s{1},1),1);
ws = zeros(size(s{1},1),1);

for j = 1:size(s{1},1)
    
    w_poly = zeros(size(gcf,2),1);
    for k = 1:size(gcf,2)
        w_poly(k) = s{1}(j,1)^(k-1);
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
    g(j,:) = (s{1}(j,1) * sh)./price';
    
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
        sh(k+1) = exp(min(cf(k,:)*w_poly,400));
    end
    sh = sh/(1+sum(sh));
    sh(1) = 1 - sum(sh);

    %get actual consumption
    cns = (egr(j) * sh)./price';

    fu = sum(bsxfun(@times,p(:,1),log(bsxfun(@plus,p(:,2),cns))),1);
    
    diff = bsxfun(@minus,cns',g);
    wp_int_int = max(normpdf(diff,zeros(size(diff)),repmat(v,size(diff,1),1)),1e-25);
    wp_int = bsxfun(@times,wp_int_int,ws)*1e7; 
    wp = prod(wp_int,2)/sum(prod(wp_int,2));
    ex = s{2}*sum(s{1}(:,2).*su.*wp); %quadrature weighted integration
    
    tu = tu + (1-alp)*fu + alp * ex;
    
    %get fundy grad
    dudc = p(:,1)./bsxfun(@plus,p(:,2),cns);
    dcds = egr(j)./price';
    sh_stack = repmat(sh(2:end)',size(cf,2),1); %once this is vectorized, it will give me the right number of repeated shares
    dsdcf = -bsxfun(@times,repmat(sh',size(cf(:),1),1).*repmat(w_poly,28,size(sh,1)),sh_stack(:)); %size(#coefs,#goods)
    for k = 1:size(cf,2)
        dsdcf(size(dsdcf,1)+k:size(dsdcf,1)+size(cf,2):end) = dsdcf(size(dsdcf,1)+k:size(dsdcf,1)+size(cf,2):end)' + sh(2:end)*w_poly(k);
    end
    dcdcf = bsxfun(@times,dcds',dsdcf);
    dudcf = dcdcf*dudc;
    grad = dudcf;
    
    %get integral grad
    sot = bsxfun(@times,wp,ones([size(s{1},1),size(s{1},1),29]));
    sot = bsxfun(@times,wp,ones([size(s{1},1),size(s{1},1),29])).*permute(sot,[2 1 3]); %this gives us wp*wp' repeated 29 times in a third direction
    sot = bsxfun(@times,permute(sot,[2 3 1]),bsxfun(@rdivide,diff,(v.^2)));
    %sot = permute(sot,[2 3 1]).*repmat(diff./repmat((v.^2),size(s{1},1),1),[1 1 size(s{1},1)]);
    sot = permute(sot,[3 1 2]);
%     fot_test = zeros(size(sot));
%     for k = 1:size(s{1},1)
%         fot_test(k,k,:) = wp(k)*-diff(k,:)./(v.^2);
%     end
%     toc
    fot = zeros(size(sot));
    fot_int = bsxfun(@rdivide,bsxfun(@times,wp,-diff),v.^2);
    ind = bsxfun(@plus,(1:size(s{1},1)+1:size(s{1},1)^2)',(0:28).*size(s{1},1)^2); 
    fot(ind) = fot_int(:);
    dwdd = fot + sot;
 
    dedw = s{2}*s{1}(:,2).*su;
     
     %expensive for loop--tried to vectorize, but it was actually slower to use arrayfun 
     dedc = zeros(29,1);
     for k = 1:29
         dedc(k) = sum((dwdd(:,:,k)'*dedw)); %simpler and equivalent way of mutliplying by dddc--and slower.... 
     end
    
    dedcf = dcdcf*dedc;
    
    grad_ex = dedcf;
    
    tg = tg + (1-alp)*grad + alp*grad_ex;
end

u = -tu;
foc = -tg;
 for k = 1:size(cf,2)
    foc(k:size(cf,2):end) = foc(k:size(cf,2):end)/scl^(k-1); 
 end
 if isnan(u) == 1
     display('NaNs off the starboard bow!  All men into lifeboats!');
 u(isnan(u) == 1 | abs(u) == inf) = 0;
 foc(isnan(foc) == 1 | abs(foc) == inf) = 0;
 end
 %scale
 u = u/pscale;
 foc = foc/pscale;
 
end