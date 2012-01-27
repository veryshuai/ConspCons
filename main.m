%temporary main file for conspicuous consumption estimation
clc;
cd '/gpfs/home/dcj138/work/ConspCons/'

%read %create data files

load dat

%open parallel functionality
matlabpool open 6

%open diary
diary 1-27-12-3.txt

%linearly extend prices to cover zeros and include only relevant price years
price = price_fixing(price); 
price = price([1,7:12,14:24],:);


%create total expenditures
exp = sum(cons,2);

c = dataset(exp,cons,char);

%drop obs if expenditure less than one dollar/yr
ind = find(exp>1);
less_than_ones = size(exp,1) - size(ind,1);
c = c(ind,:);

%drop obs if any good has negative expenditure
neg_exps = 0;
for k = 1:29
ind = find(c.cons(:,k)>=0);
c = c(ind,:);
neg_exps = neg_exps + size(c.cons(:,k),1)-size(ind,1);
end
display(['Lost ', num2str(less_than_ones), ' expenditures less than one dollar, and further lost ', num2str(neg_exps),' negative expenditure observations']); 

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

% create expenditure grid by year [year,grid]
gr_siz = 10;
w_b = zeros(2,18);
egr = zeros(18,gr_siz);
for k = 1:18
w_b(:,k) = quantile(c.exp(c.y==k),[.1;.9]);
egr(k,:) = logspace(log(w_b(1,k))/log(10),log(w_b(2,k))/log(10),gr_siz);
end

%create type indices
ti = zeros(size(c,1),4);
ti(:,1) = c.char(:,3)==1 & c.char(:,4) == 1; %bo4
ti(:,2) = c.char(:,3)==1 & c.char(:,4) == 0; %no4
ti(:,3) = c.char(:,3)==0 & c.char(:,4) == 1; %bu4
ti(:,4) = c.char(:,3)==0 & c.char(:,4) == 0; %nu4

%make a type variable
c.type = ti(:,1)+2*ti(:,2) +3*ti(:,3)+4*ti(:,4); 

%create two dimensional consumption array with dimensions year, type, and element expenditure level
c = sortrows(c,{'exp'});
ca = cell(18,4);
ce = cell(18,4);
for k = 1:18
    for m = 1:4
       ind = find(c.y==k & c.type==m);
        ca{k,m} = c.cons(ind,:);
        ce{k,m} = c.exp(ind);
    end
end

% find the cutoff observations in data, this is going to be a two dimensional array [year,type]
% and the elements will be the cuts 
egr_cuts = cell(18,4);
for k = 1:18
    for m = 1:4
        egr_cuts{k,m} = zeros(size(egr,2),1);
        for j=1:size(egr,2)
            egr_cuts{k,m}(j) = find(ce{k,m}>egr(k,j),1,'first');
        end
    end
end

%wealth dist calculations, matrix [density,year,type]
w = zeros(size(egr,2),18,4);
for j = 1:18
    for k = 1:4  
        w(:,j,k) = histc(ce{j,k},egr(j,:))+1e-6; %added a bit to ensure that there are no zeros (full support)
        w(:,j,k) = w(:,j,k)/sum(w(:,j,k),1); %make it a dist
    end
end

%v is literally the inverse of the vindex times grid size
v = zeros(4,29);
for k = 1:4
    v(k,:) = vin(k,:).^-1;
end

%call genetic algorithm

pop = [];

options=gaoptimset('Display','iter','PopulationSize',18,'Generations',150,... 
   'StallTimeLimit',86400,'TimeLimit',Inf,'MutationFcn',@mutationadaptfeasible,...
   'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
   'PlotFcns',@gaplotbestf,'EliteCount',1);

    [X,fval,exitflag,output,population,scores] = ga(@(X) likelihood(X,price,egr,v,w,ca,egr_cuts),60,...
    [],[],[],[],[ones(1,29)*1e-12,ones(1,29)*1,0,1e2],[ones(1,29)*2,ones(1,29)*2000,1,1e4],[],options);  

matlabpool close
diary off



