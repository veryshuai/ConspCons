%temporary main file for conspicuous consumption estimation
clc;
cd '/gpfs/home/dcj138/work/ConspCons/'

%read %create data files

load dat

price = price_fixing(price); %linearly extend prices to cover zeros

%create total expenditures
exp = sum(cons,2);

c = dataset(cons,char,exp);
%drop obs if expenditure less than one dollar/yr
ind = find(exp>0);
c = c(ind,:);
%drop obs if any good has negative expenditure
for k = 1:29
ind = find(c.cons(:,k)>=0);
c = c(ind,:);
end

egr = linspace(1,max(c.exp),1e4);

%for now I pool years and do the wealth dist calculation
w = histc(c.exp,egr)+1e-6;
w = w/sum(w); %make it a dist

%v is literally the inverse of the vindex, using whites under40
v = vin(1,:).^-1;

%as an initial value, use the simple stone geary demand system
param1 = ones(29,1);
param2 = ones(29,1);
sgd = bsxfun(@plus,-param2',bsxfun(@times,param1'./price(1,:),egr'+sum(param2.*price(1,:)'))/sum(param2));
sgd(sgd<0) = 0; %corner

%for now let the unconstrained thang be the social guess.
g = sgd;

%now lets add optimization given prices.  to do this I need to figure out
%the first order condions given g.

options=optimset('Display','iter','jacobian','on','TolFun',1e-6,'TolX',1e-8,'DerivativeCheck','off',...
    'GradObj','on');

c_mat = zeros(size(g));
c_mat(1,:) = ktrlink(@(x) util_FOC(ones(29,2),x,g,.5,v,w,1),log(ones(29,1)/29),[],[],ones(1,29),egr(1),[],[],[],options);
%c_mat(1,:) = ktrlink(@(x) -util_FOC(ones(29,2),x,g,.5,v,w,1),[log(ones(29,1)/29);1],[],[],[],[],[],[],[],options);
for k = 2:size(egr,2)
    c_mat(k,:) = ktrlink(@(x) util_FOC(ones(29,2),x,g,.5,v,w,egr(k)),c_mat(k-1,:)',[],[],ones(1,29),egr(1),[],[],[],options);
end
c_mat = exp(c_mat);




