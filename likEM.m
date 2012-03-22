% This calls the likelihood in the new EM version of my ConspCons project

clear;

cd '/gpfs/home/dcj138/work/ConspCons/'

%read %create data files

load dat

%open parallel functionality
%if matlabpool('size')<12
%matlabpool open 12
%end

%open diary (NOW USING LOG FILES called in matlab command line when opening)
%diary 2-28-12.txt

%linearly extend prices to cover zeros and include only relevant price years
price = price_fixing(price); 
price = price([1,7:12,14:24],:);

%create total expenditures
exp = sum(cons,2);

c = dataset(exp,cons,char);

%drop obs if expenditure less than one dollar/yr
ind = find(exp>1000);
less_than_ones = size(exp,1) - size(ind,1);
c = c(ind,:);

%drop obs if any good has negative expenditure
neg_exps = 0;
for k = 1:29
ind = find(c.cons(:,k)>=0);
c = c(ind,:);
neg_exps = neg_exps + size(c.cons(:,k),1)-size(ind,1);
end
display(['Lost ', num2str(less_than_ones), ' expenditures less than one thousand dollars, and further lost ', num2str(neg_exps),' negative expenditure observations']); 

%create year numbers
c.y = zeros(size(c.exp));
c.y(c.char(:,7)==80) = 1;
for k = 2:7
    c.y(c.char(:,7)==84+k) = k;
end
for k = 8:14
    c.y(c.char(:,7)==85+k) = k;
end
for k = 15:18
    c.y(c.char(:,7)==k-15) = k;
end

%Get expenditure lower bounds by year 
w_b = zeros(18,1);
c_b = zeros(18,1);
for k = 1:18
w_b(k) = min(c.exp(c.y==k));
c_b(k) = max(max(c.cons(c.y==k)));
end

%create type indices
ti = zeros(size(c,1),4);
ti(:,1) = c.char(:,3)==1 & c.char(:,4) == 1; %bo4
ti(:,2) = c.char(:,3)==1 & c.char(:,4) == 0; %no4
ti(:,3) = c.char(:,3)==0 & c.char(:,4) == 1; %bu4
ti(:,4) = c.char(:,3)==0 & c.char(:,4) == 0; %nu4

%make a type variable
c.type = ti(:,1)+2*ti(:,2) +3*ti(:,3)+4*ti(:,4); 

%vindex probabilities
v = cell(4,1);
for k = 1:4
    v{k} = vin(k,:)/sum(vin(k,:));
end

%make random observation type index
%ot_ind = randi(29,size(c,1),1);
%ot_ind_old = ones(size(ot_ind));

%make c of a managable size
%I = randi(size(c,1),10000,1);
load I;
c = c(I,:);
%ot_ind = ot_ind(I);
load ot_ind;
ot_ind_old = ot_ind_old(I);

%change to regular arrays (seems to be faster than dataset)
cons = c.cons;
exp = c.exp;
y = c.y;
type = c.type;

init = [1;0.888500119023610;0.938173068223820;0.883602537604718;0.950723378567813;1.07876791156688;0.606974395481645;0.918292576001753;0.936254014486101;1.09078531261364;0.866613598387630;1.08811339132943;0.863571200630444;0.881897486236716;1.34105044786735;0.876096616765889;1.28375646624514;1.98072496565742;0.914572775181359;0.952244899959159;0.564885647510619;0.912772934257350;0.863122227812932;0.897935087167247;0.896084706640516;0.711819865715615;1.02752969128126;0.892146822623657;0.950221981050088;6.69470868590175;6.68840200751367;6.68477285991830;6.68857507346183;6.68834383646738;6.70022827930813;6.59586708101165;6.68848160687392;6.68835625633765;6.68468460494704;6.68467121128477;6.68952940924100;6.68836091640282;6.68834891580383;6.68799193303063;6.68831082517213;6.68846772879441;6.68928783323428;6.68811109401570;6.68804634665044;6.71043625473747;6.68833796230337;6.68835017385212;6.68833378011546;6.68842671576891;6.71062769592365;6.68793973621889;6.68839374434560;6.68840856073787;0.0561122405210935;];

%% EM algorithm
tic
while min(ot_ind == ot_ind_old) == 0

% E step

options=optimset('Display','iter','UseParallel','always','Algorithm','interior-point','TolCon',1e-3,'TolX',1e-6,'TolFun',1e-3,'MaxFunEvals',1e5,'MaxIter',10,'SubProblemAlgorithm','cg');%,'Hessian','bfgs');

[X,fval,exitflag] = ktrlink(@(x) likEM_loop(x,cons,exp,y,price,ot_ind,w_b,c_b),init,[],[],[],[],...
    [1;zeros(28,1);ones(29,1)*0;0],[1;ones(28,1)*2;ones(29,1)*inf;1],...
    [],options);   

init = X;
display(X);

% M step

ot_ind_old = ot_ind;
ot_ind = not_ind(init,cons,exp,y,type,price,w_b,c_b,v);

if toc>2*3600
   break; 
end

end

display('All Done!');

save 3-21-12-1.mat


