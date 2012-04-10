% This script creates plots for use in the writeup of my conspicuous consumption
% paper 

%% Preliminary stuff

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

%% Which year?
% how about 2001?
I = find(c.y==16);

%% Create engels curve plots

colors = vertcat(repmat([1,.5,.5],size(1:round(9*size(I,1)/10),2),1),repmat([1,.25,.25],size(round(9*size(I,1)/10)+1:round(39*size(I,1)/40),2),1),repmat([1,0,0],size(round(39*size(I,1)/40)+1:size(I,1),2),1));

shr = figure(1);
%title('Log Expenditure Shares by Log Expenditure');
num = 0;
lab = {'food at home','food out','tobacco','alcohol at home','alcohol out','clothing','laundry','jewelry','rent','home','hotel','furniture','utilities','telephone','home insurance','medical care'};
for k = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
   num = num+1;
   subplot(4,4,num) 
   scatter(log(c.exp(I)),log(c.cons(I,k)./c.exp(I)),2,colors,'filled');
   xlabel(lab{num})
   axis([6 14 -15 0]);
   set(gca,'FontSize',8);
end
print(shr,'-dpdf','figures/shares');

lev = figure(2);
%title('Log Expenditure Shares by Log Expenditure');
num = 0;
lab = {'food at home','food out','tobacco','alcohol at home','alcohol out','clothing','laundry','jewelry','rent','home','hotel','furniture','utilities','telephone','home insurance','medical care'};
for k = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
   num = num+1;
   subplot(4,4,num) 
   scatter(log(c.exp(I)),log(c.cons(I,k)),2,colors,'filled');
   xlabel(lab{num})
   axis([6 14 0 11]);
   set(gca,'FontSize',8);
end
print(lev,'-dpdf','figures/levels');

%% Create Wealth Histogram

exphist = figure(3);

hist(log(c.exp(I)),100);
xlabel('Log Expenditure');
ylabel('Frequency');
title('Histogram of Log Expenditure over Consumers, 2001'); 
print(exphist,'-dpdf','figures/exphist');