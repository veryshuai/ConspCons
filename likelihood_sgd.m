function r = likelihood_sgd(X,price,egr,v,w,ca,egr_cuts)

warning off all;

%tic

param1 = X(1:29)'; %cobb douglas parameters
param2 = X(30:58)'; %stone geary parameters
%alp = X(59); %weight of conspicuous consumption
%vm = X(60); %v multiplier

%alp = 0;

%v = v*vm; %scale v by the multiplier

%as an initial value, use the simple stone geary demand system

options=optimset('Display','off','jacobian','on','TolFun',1e-3,'TolX',1e-5,'TolCon',1e-6,'DerivativeCheck','off',...
    'GradObj','on');%,'Algorithm','active-set');

%get "naive" initial value
%init = zeros(size(egr,2),29,size(egr,1));
%for k = 1:18
%init(:,:,k) = bsxfun(@plus,-param2,bsxfun(@times,param1/sum(param1)./price(k,:)',(egr(k,:)'+sum(param2.*price(k,:)'))'))';
%end
%init(init<=0) = 0; %corner
%for k = 1:18
%    init(:,:,k) = bsxfun(@times,egr(k,:)',bsxfun(@rdivide,init(:,:,k),sum(bsxfun(@times,price(k,:),init(:,:,k)),2)));
%end
    
tic 

pt = 3; %number of polynomial terms

scl = max(egr(:));
init = [-0.876908208149109;-0.388138631036905;0.489678136868298;-2.17661927161956;0.973990796532659;1.04852115604715;-5.40323995158061;0.330086137471565;0.910529521672106;-3.33257241473697;0.871196554392692;1.06816131627565;-3.70844073042796;0.816392447513688;1.05814597945497;-1.41761141677479;-0.150470463449478;0.599010852085380;-4.04977337269676;0.892820478102040;1.10186759972587;-6.61572556806999;0.0603469476663820;0.837695845464407;-4.97281548403554;0.704569132131087;1.05346332285584;0.222054621323285;-1.63456919605491;0.225842675862680;-0.841754910320384;-0.782423097832529;0.364518426320688;-0.665442785573945;-0.900628188761481;0.319999663050599;-9.20201927107749;-0.721246600123931;0.594836555637205;-0.698528904329736;-1.52297162202306;0.135052345902352;-1.03862944162559;0.295278224069180;0.722809168264034;-1.34436128705644;0.0659961337259723;0.668735756732909;-1.85886927622982;0.855996831182266;0.983844346698051;-8.42927606519891;-0.462841617575564;0.678002475435340;-1.56649366702422;2.57248242044950;1.58143060916988;-2.66596538362298;1.12350904592443;1.13504704831237;-3.01223686886268;0.908611954069789;1.06858968779075;-1.71925032300938;1.00526401211105;1.02816906557097;-7.43526865558882;-0.312544134884067;0.708419444535324;-3.04892092088874;1.16770048520396;1.18112339106101;-9.61639032377301;-0.856541254826385;0.551860833219124;-1.75165301012071;1.82985509667552;1.33948334567094;-1.01263451856510;0.431058131347962;0.764699901238742;-9.48110443690475;-0.818895309618202;0.563125596313716;];

sgd = zeros(size(init,1)/28,28,size(egr,1));
for m = 1:18
    temp_resize = fminunc(@(x) util_FOC_sgd([param1,param2],x,egr(m,:,1),price(m,:),scl),init,options);
    init = temp_resize;
    sgd(:,:,m) = reshape(temp_resize,size(init,1)/28,28);
end

%implied expenditures for each wealth type
expend = zeros(size(egr,2),29,18);
w_poly = zeros(pt,1);
share = zeros(size(egr,2),29,18);
sh = zeros(29,1);
for m = 1:18
    for j = 1:size(egr,2)
        for k = 1:pt
            w_poly(k) = (egr(m,j,1)/scl)^(k-1);
        end
        for k = 2:29
            sh(k) = exp(sgd(:,k-1,m)'*w_poly);
        end
        sh = sh/(1+sum(sh));
        sh(1) = 1 - sum(sh);
        share(j,:,m) = sh;
        expend(j,:,m) = egr(m,j,1)*sh./price(m,:)'; 
    end
end

%put in types for resid calculations
g = zeros(size(expend,1),size(expend,2),size(expend,3),4);
for k = 1:4
    g(:,:,:,k) = expend;
end

% init = g;
% c_mat = zeros(size(init));
% 
% options = optimset('Display','iter');
% 
% for m = 1:18
%    for j = 1:4
%        for l = 1:size(egr,2)
%         [c_mat(l,:,m,j),~,flag] = fminsearch(@(x) inner(x,[param1,param2],alp,v,w,price,egr),init,options); 
%        end
%    end
% end

%INNER LOOP
%find a fixed point such that actual consumption matches expectataions

% options=optimset('Display','off','jacobian','on','TolFun',1e-4,'TolX',1e-4,'DerivativeCheck','off',...
%     'GradObj','on','Algorithm','active-set');
% 
% err = 1;
% err_lag = inf;
% g_lag = g;
% count = 0;
% punishment = 0;
% break_flag = 0;
% while err > 1e-1 && count<20
%     count = count+1;
%     %display(count);
%     l_count = 0;
%     c_mat = zeros(size(g));
%     for m = 1:18    
%             l_count = l_count+1;
%             %display(l_count)
%             for j = 1:4
%                 for k = 1
%                     [c_mat(k,:,m,j),~,flag] = ktrlink(@(x) util_FOC([param1,param2],x,g(:,:,m,j),alp,v(j,:),w(:,m,j)),g(k,:,m,j),[],[],price(m,:),egr(m,k),zeros(29,1),[],[],options);
%                     if flag ~= 0 && flag ~= -100
%                     %[c_mat(k,:,m,j),~,flag] = fmincon(@(x) util_FOC([param1,param2],x,g(:,:,m,j),alp,v(j,:),w(:,m,j)),g(k,:,m,j),[],[],price(m,:),egr(m,k),zeros(29,1),[],[],options);
%                     %if flag < 1
%                         %display(flag);
%                         punishment = 1;
%                     end
%                 end
%                 for k = 2:size(egr,2)
%                     [c_mat(k,:,m,j),~,flag] = ktrlink(@(x) util_FOC([param1,param2],x,g(:,:,m,j),alp,v(j,:),w(:,m,j)),g(k,:,m,j),[],[],price(m,:),egr(m,k),c_mat(k-1,:,m,j),[],[],options);%zeros(29,1),[],[],options);
%                     if flag ~= 0 && flag ~= -100
%                     %[c_mat(k,:,m,j),~,flag] = fmincon(@(x) util_FOC([param1,param2],x,g(:,:,m,j),alp,v(j,:),w(:,m,j)),g(k,:,m,j),[],[],price(m,:),egr(m,k),zeros(29,1),[],[],options);%c_mat(k-1,:,m,j),[],[],options);
%                     %if flag < 1
%                         %display(flag);
%                         punishment = 1;
%                     end
%                     if isnan(c_mat(k,1,m,j)) == 1
%                         display('uh oh, we have NaNs');
%                         break_flag = 1;
%                         punishment = 1;
%                         break;
%                     end
%                 end
%                 if break_flag == 1
%                     display('NaN link 1');
%                     break;
%                 end
%             end
%             
%         if break_flag == 1
%             display('NaN link 2');
%             break;
%         end
%     end
%     err = norm(g(:)-c_mat(:));
%     if err <= err_lag
%         g_lag = g;
%         g = c_mat;
%         err_lag = err;
%     else
%         g = (g + g_lag)/2;
%     end
%     %display(err_lag);
%     %display(err);
%     if count == 20
%             punishment = 1;
%     end
%     if break_flag == 1  || isnan(err) == 1
%         display('NaN link 3');
%         c_mat = zeros(size(g));
%         break;
%     end    
% end

display(param1);
display(param2);
%display(alp);
%display(vm);

% c_mat_vals = zeros(size(g));
% for k = 1:size(egr,2)
%     for m = 1:18
%         for j = 1:4
%             c_mat_vals(k,:,m,j) = g(k,:,m,j).*price(m,:);
%         end
%     end
% end

r = -resid(g,ca,egr_cuts);
%if punishment == 1 
%    display('Punishment added for convergence issue.');
%    r = r*1.5; %add a 50% punishment for failure to converge
%end
display(r);

%g = zeros(size(sgd,1),size(sgd,2),size(sgd,3),4);
%for k = 1:4
%    g(:,:,:,k) = sgd;
%end

%r_test = -resid(g,ca,egr_cuts);
%display(r_test);

%adv = r_test-r; %advantage of including consp cons in terms of fit
%max_diff = display(max(max(max(max(abs(g-c_mat))))));  %maximum difference between consumption over regime type

%display(adv);
%display(max_diff);

%toc

%diary off
diary 2-4-12-sgd-1.txt

end