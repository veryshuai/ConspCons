%This script plots histograms of the consumption data.

%clear;

cd '/gpfs/home/dcj138/work/ConspCons/'

%read %create data files

load dat

share = bsxfun(@rdivide,cons,sum(cons,2));
lab = {'food at home','food out','tobacco','alcohol at home','alcohol out','clothing','laundry','jewelry','rent','home','hotel','furniture','utilities','telephone','home insurance','medical care'};

shrplot = figure(1);
for k = 1:16
    h = subplot(4,4,k);
    hist(log(share(share(:,k)>0,k)),75);
    %set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'FontSize',5);
    xlim([-15 5]);
    xlabel(lab{k})
end
print(shrplot,'-dpdf','figures/shrplot');

zeromass = zeros(29,1); 

for k = 1:29
   zeromass(k) = sum(cons(:,k) == 0)/size(cons,1);
end

display(zeromass);