% This calls the likelihood in the new EM version of my ConspCons project

%clear;

cd '/gpfs/home/dcj138/work/ConspCons/'

%read %create data files

load dat

%open parallel functionality
%if matlabpool('size')<12
%    matlabpool open 12
%end

%open diary (NOW USING LOG FILES called in matlab command line when opening)
diary 4-7-12.txt

%linearly extend prices to cover zeros and include only relevant price years
price = price_fixing(price); 
price = price([1,7:12,14:24],:);

%create total expenditures
expend = sum(cons,2);

c = dataset(expend,cons,char);

%drop obs if expenditure less than one dollar/yr
ind = find(expend>1000);
less_than_ones = size(expend,1) - size(ind,1);
c = c(ind,:);

%drop obs if any good has negative expenditure
neg_exps = 0;
for k = 1:29
ind = find(c.cons(:,k)>=0);
neg_exps = neg_exps + size(c.cons(:,k),1)-size(ind,1);
c = c(ind,:);
end
display(['Lost ', num2str(less_than_ones), ' expenditures less than one thousand dollars, and further lost ', num2str(neg_exps),' negative expenditure observations']); 

%drop obs if food at home is zero
ind = find(c.cons(:,1)>0);
nofood = size(c.cons,1)-size(ind,1);
c = c(ind,:);
display(['Also lost ', num2str(nofood), ' because of zero expenditures on food at home']); 
display(['We have a total of ', num2str(size(c.expend,1)), ' usable observations']);

%create year numbers
c.y = zeros(size(c.expend));
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
w_b(k) = min(c.expend(c.y==k));
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
ot_ind = randi(29,size(c,1),1);
ot_ind_old = ones(size(c,1),1);

%make c of a managable size
I = randi(size(c,1),5000,1);
%load I;
c = c(I,:);
ot_ind = ot_ind(I);
%load ot_ind;
ot_ind_old = ot_ind_old(I);

%change to regular arrays (seems to be faster than dataset)
cons = c.cons;
expend = c.expend;
y = c.y;
type = c.type;

%init = [1;7.59210648140400;4.73543180019105;5.74862875334654;4.95440784574850;12.7037395527745;5.69478318955006;7.88753934896234;1.44938461037568;9.18538421369885;9.05132228657440;9.00566247190964;3.49072874128893;8.27878887702105;10.4018645159034;8.17370675069044;15.4762574500119;9.68697383660247;3.55478951284719;7.94333412618397;7.18239264327919;11.4795764405883;4.45354347075682;9.18123307895688;4.89033935794127;12.1057696447419;9.14687789130620;9.31117750145790;7.41878093944565;7.61952386298256;7.70619388603037;8.23975939816500;7.71149476704346;7.84473732613889;7.25724006498097;7.48925368453039;7.50087264872762;7.82344615001618;7.60552520095296;7.57603338714131;7.19619527149678;7.87225119947359;7.39836937524840;7.30019148841177;7.81848287278297;7.73258241738970;9.12571047476629;7.32590778396217;7.73716194637371;8.15061049349470;7.69021491224891;7.67405438242449;7.66479184955140;7.72149575366749;6.55816617452169;7.71298023499229;7.98213121498031;7.23763544816502;0.438919082280316;];
%init = linspace(0,1,36)';
%init = [.5,ones(1,56)];
init = .0477;
% init = [
%      %0.1 0.1714    1.1580    1.6976    0.1727   -0.3817   -0.7895    1.7836  0.0798    1.7731    1.8077    1.6738    0.2052    0.8772    0.7889 1.2605    0.3888    1.3579    1.2662   -0.3652   -0.2960    0.0327 0.1665 -0.3678    0.7264    0.3889    0.2509    0.6431 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;...
%      0.451 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.401 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.351 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.251 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.151 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.051 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.101 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.5 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      0.0 -0.0786    0.8455    0.1976    0.2821   -1.6317   -0.7895    0.2505  -1.6702    1.0231    0.7239    0.7344   -0.5448   -1.1228    1.0389 0.0727   -0.3783    1.8735   -0.5010   -1.3339   -0.2960   -0.9673 0.1665   -1.3678    1.1014    0.3889   -1.2491   -0.0444    0.4552 2.0000    3.5000    4.5000    3.9375    1.8438    3.1250    4.0937 2.2500    2.0781    4.2500    3.2969    1.5000    1.5000    4.2500 2.6179    5.1094    5.0625    3.2500    1.6605    2.0156    0.8585 4.3750    3.0625    5.3125    3.5000    1.3532    3.5000    3.8750 
%      ];
%scale
%scl = abs(likEM_loop(init,cons,expend,y,price,ot_ind,w_b,c_b,1));
scl = 1;

%Calculate zero mass point probabilities before hand
pc = cons>0;
zm = sum(pc)/size(pc,1);
zm = 1-zm(2:end);

%Get initial ot_ind
ot_ind = not_ind_oa(.05*ones(57,1),cons,expend,y,type,price,ot_ind,w_b,c_b,scl,v,zm');
%% EM algorithm
tic
totiters = 0;
while min(ot_ind == ot_ind_old) == 0

distance = norm(ot_ind-ot_ind_old)/norm(ot_ind_old);    
totiters = totiters+1;
display(totiters);
display(distance);

% E step

%options=optimset('Display','iter','Algorithm','interior-point','TolCon',1e-10,'TolX',1e-12,'TolFun',1e-10,'MaxFunEvals',1e5,'MaxIter',50);%,'SubProblemAlgorithm','cg');%,'Hessian','bfgs');
%options=optimset('Display','iter','LargeScale','off','TolCon',1e-12,'TolX',1e-12,'TolFun',1e-10,'MaxFunEvals',1e5,'MaxIter',150,'HessUpdate','bfgs');%,'SubProblemAlgorithm','cg');%,'Hessian','bfgs');
options = gaoptimset('Display','iter','InitialPopulation',init,'PopulationSize',12,'UseParallel','always','PlotFcns',@gaplotbestf,'StallGenLimit',10);

%[X,fval,exitflag] = fminunc(@(x) likEM_loop(x,cons,expend,y,price,ot_ind,w_b,c_b,scl,v),init,options);
X = ga(@(x) likEM_loop_oa(x,cons,expend,y,price,ot_ind,w_b,c_b,scl,v,type,zm),1,[],[],[],[],[0,-2*ones(1,28),zeros(1,28)],[.5,3*ones(1,28),6*ones(1,28)],[],[],options);   
%X = ktrlink(@(x) likEM_loop(x,cons,expend,y,price,ot_ind,w_b,c_b,scl,v,type,zm),init,[],[],[],[],[0,-inf*ones(1,28),zeros(1,28)],[1,inf*ones(1,56)],[],options);   

[~,mu,sig,zm] = likEM_loop_oa(X,cons,expend,y,price,ot_ind,w_b,c_b,scl,v,type,zm);

init = X;
%implied_alp = exp(X)/(1+exp(X));
implied_alp = X;

display([mu,sig,zm]);
display(implied_alp);

% M step

ot_ind_old = ot_ind;
ot_ind = not_ind_oa([X;mu;sig],cons,expend,y,type,price,ot_ind,w_b,c_b,scl,v,zm);

%if toc>2*3600
%   break; 
%end

end

display('All Done!');

save 4-7-12.mat


