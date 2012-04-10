% This calls the likelihood in the new EM version of my ConspCons project

clear;

cd '/gpfs/home/dcj138/work/ConspCons/'

%read %create data files

load dat

%create total expenditures
expend = sum(cons,2);

wp = lognfit(expend(expend>0));

load results

%open parallel functionality
%if matlabpool('size')<12
%    matlabpool open 12
%end

%open diary (NOW USING LOG FILES called in matlab command line when opening)
%diary 3-31-12-2.txt

%linearly extend prices to cover zeros and include only relevant price years
price = price_fixing(price); 
price = price([1,7:12,14:24],:);

%vindex probabilities
v = cell(4,1);
for k = 1:4
    v{k} = vin(k,:)/sum(vin(k,:));
end
cumv = v{2};
cumv = cumv(2:end)/sum(cumv(2:end));
cumv = cumsum(cumv);

%simulate individuals
  siz = 5000;
%  drw = rand(siz,28,2);
%  wdrw = rand(siz,1);
%  otdrw = rand(siz,1);
%  clearvars -except drw wdrw otdrw
%  save shocks/drw;

load shocks/drw;

% for now use this:
%mu = ...
%    [0.3509    0.4813    0.7607    0.8393    0.1570    0.9851    0.5268 0.6842    0.5361    0.5660    0.8003    0.6397    0.3393    0.2741 0.7511    0.6939    0.6105    0.0312    0.8866    0.1638    0.2360 0.2003    0.2593    0.2902    0.6449    0.3176    0.0595    0.4525];

%sig = ...
%    [0.5971    0.6922    0.7841    0.9373    0.9684    0.1811    0.2434 0.2394    0.8141    0.2226    0.7321    0.4354    0.3959    0.6472 0.1486    0.5036    0.6386    0.8005    0.7983    0.4237    0.4660 0.2926    0.5449    0.7544    0.6317    0.8625    0.6405    0.3878];

%alp = .1610;

%zm = sum(cons>0,1)./sum(cons>=0,1);
%zm = zm(2:end);

mp = drw(:,:,1)>repmat(zm',siz,1);
gam = mp.*logninv(drw(:,:,2),repmat(mu',siz,1),repmat(sig',siz,1));
gam = horzcat(ones(size(gam,1),1),gam);
w = logninv(wdrw,wp(1),wp(2));
ot = zeros(siz,1);
for k = 1:siz
    ot(k) = 1 + find(otdrw(k)<cumv,1,'first');
end

%find utility maximizing taxes

init = ones(29,1)*1.01;

options = optimset('Display','iter');

tax = fmincon(@(x) ag_util(x,price(18,:),gam,w,ot,alp),init,[],[],[],[],ones(29,1),ones(29,1)*3,[],options);
save aftertax
[~,util_new,adjc_new] = ag_util(tax,price(18,:),gam,w,ot,alp);
[~,util_old,adjc_old] = ag_util(ones(29,1),price(18,:),gam,w,ot,alp);
welfgain = (util_new(util_new>0)-util_old(util_new>0))./util_old(util_new>0);
welfgain(abs(welfgain)>.2) = 0;
labs = ones(size(ot(util_new>0)))*NaN;
I = randi(size(ot(util_new>0),1),0,1);
labs(I) = 1;
%labs = ones(size(ot(util_new>0)));
taxscat = figure(3);
lscatter(log(w(util_new>0)),welfgain,ot(util_new>0).*labs,'MissingLabel','*','TextColor','black','FontSize',10);
print(taxscat,'-dpdf','figures/taxscat');

%% Plot the "fake" data just as I did the real data, but in BLUE

adjw_old = sum(adjc_old,2);

colors = vertcat(repmat([.5,.5,1],size(1:round(9*size(adjw_old,1)/10),2),1),repmat([.25,.25,1],size(round(9*size(adjw_old,1)/10)+1:round(39*size(adjw_old,1)/40),2),1),repmat([0,0,1],size(round(39*size(adjw_old,1)/40)+1:size(adjw_old,1),2),1));

shr = figure(4);
%title('Log Expenditure Shares by Log Expenditure');
num = 0;
lab = {'food at home','food out','tobacco','alcohol at home','alcohol out','clothing','laundry','jewelry','rent','home','hotel','furniture','utilities','telephone','home insurance','medical care'};
for k = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
   num = num+1;
   subplot(4,4,num) 
   scatter(log(adjw_old),log(adjc_old(:,k)./adjw_old),2,colors,'filled');
   xlabel(lab{num})
   axis([6 14 -15 0]);
   set(gca,'FontSize',8);
end
print(shr,'-dpdf','figures/shares_fake');

lev = figure(5);
%title('Log Expenditure Shares by Log Expenditure');
num = 0;
lab = {'food at home','food out','tobacco','alcohol at home','alcohol out','clothing','laundry','jewelry','rent','home','hotel','furniture','utilities','telephone','home insurance','medical care'};
for k = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
   num = num+1;
   subplot(4,4,num) 
   scatter(log(adjw_old),log(adjc_old(:,k)),2,colors,'filled');
   xlabel(lab{num})
   axis([6 14 0 11]);
   set(gca,'FontSize',8);
end
print(lev,'-dpdf','figures/levels_fake');

%% now figure out the gains for an "average" person of each observation type from
%eliminating the motive for conspicuous consumption.

avp = exp(mu-sig.^2);

sig_util = zeros(1000,29);
nosig_util = zeros(1000,29);
wealth = linspace(1000,20000,1000)';
for k = 2:29
   [sig_util(:,k),nosig_util(:,k)] = util_calc(wealth,alp,price(18,:),avp,k);
end

dif = (nosig_util-sig_util)./sig_util;
nohouse = [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29];
uchnge = figure(1);
plot(dif(:,nohouse));
figure(2);
plot(dif(:,10));
print(uchnge,'-dpdf','figures/uchnge');

save opttax