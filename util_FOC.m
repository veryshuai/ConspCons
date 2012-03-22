function [u,foc] = util_FOC(p,cf,price,scl,gcf,alp,v,egr,pscale,P)
%calculates utility FOC for parameter vector p, approx coefficients cf and
%weath level w

%get social guess of utility for large number of wealth types (1000) and 29 observation types.
gcf = reshape(gcf,size(gcf,1)/29,29)';
for k = 1:size(gcf,2)
    gcf(:,k) = gcf(:,k)/scl^(k-1);
end

%wealth grid size and actual wealth grid
ws = 1000;
wg = linspace(egr(1),egr(end),ws);
nwg = wg-sum(price'.*p(:,2));

%This will hold a matrix of consumption guesses for observation types, each cell element is 1000x29
g  = cell(29,1);

%weighting vector
wv = p(:,1)./price'; 

%wealth powers
    w_poly = ones(ws,P);
    for k = 2:P
        w_poly(:,k) =wg.^(k-1);
    end

gc = cell(29,1);
%Loop over announced observation types
for j = 1:29
    
        %get consumptions 
	gc{j} = zeros(29,ws);
	gc{j}(j,:) = p(j,2) + (wg' - sum(price'.*p(:,2))).*(exp(sum(bsxfun(@times,gcf(j,:),w_poly),2))./(1+exp(sum(bsxfun(@times,gcf(j,:),w_poly),2))))/price(j);
	nnwg = nwg - gc{j}(j,:);
	ctemp = p(:,2) + wv*nnwg;
	ctemp(j,:) = gc{j}(j,:);

	gc{j} = ctemp;

	%get utilities
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
    wp_int = bsxfun(@times,wp_int_int,ws); 
    wp = prod(wp_int,2);%prod(wp_int,2)/sum(prod(wp_int,2));
    ex = s{2}*sum(s{1}(:,2).*su.*wp);%/sum(s{1}(:,2)); %quadrature weighted integration
    
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
 
    dedw = s{1}(:,2).*su/sum(s{1}(:,2));
     
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
