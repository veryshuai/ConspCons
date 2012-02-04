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

options=optimset('Display','off','jacobian','on','TolFun',1e-6,'TolX',1e-6,'TolCon',1e-6,'DerivativeCheck','off',...
    'GradObj','on');%,'Algorithm','active-set');

%get "naive" initial value
init = zeros(size(egr,2),29,size(egr,1));
for k = 1:18
init(:,:,k) = bsxfun(@plus,-param2,bsxfun(@times,param1/sum(param1)./price(k,:)',(egr(k,:)'+sum(param2.*price(k,:)'))'))';
end
init(init<=0) = 0; %corner
for k = 1:18
    init(:,:,k) = bsxfun(@times,egr(k,:)',bsxfun(@rdivide,init(:,:,k),sum(bsxfun(@times,price(k,:),init(:,:,k)),2)));
end
    
tic 

sgd = zeros(size(egr,2),29,size(egr,1));
for j = 1:size(egr,2)
    for m = 1:18
        [sgd(j,:,m),~] = fmincon(@(x) util_FOC_sgd([param1,param2],x,init(j,:,k),0,v(1,:),w(:,m,1)),init(k,:,m),[],[],price(m,:),egr(m,k),zeros(29,1),[],[],options);
    end
end

toc 

%let the initial guess be the sgd.
g = zeros(size(sgd,1),size(sgd,2),size(sgd,3),4);
for k = 1:4
    g(:,:,:,k) = sgd;
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

c_mat_vals = zeros(size(g));
for k = 1:size(egr,2)
    for m = 1:18
        for j = 1:4
            c_mat_vals(k,:,m,j) = g(k,:,m,j).*price(m,:);
        end
    end
end

r = -resid(c_mat_vals,ca,egr_cuts);
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
%diary on

end