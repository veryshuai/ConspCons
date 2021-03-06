%Read in Consp Cons data

clear;
clc;

cd '/gpfs/home/dcj138/work/CEX_code/'

consdat = dlmread('merged_cex.csv',',',1,0);

mainspread = dlmread('MainSpread.csv',',',1,2);

%Categorizing all the things:

%cons: each row is an observation, and expenditures in the 29 categories 
%cons order: 1.FdH, 2.FdO, 3.Cig, 4.AlH, 5.AlO, 6.Clo, 7.Lry, 8.Jwl, 9.Brb,
%10.Hom, 11.Htl, 12.Fur, 13.Utl, 14.Tel, 15.HIn, 16.Med, 17.Fee, 18.LIn, 
%19.Car, 20.CMn, 21.Gas, 22.CIn, 23.Bus, 24.Air, 25.Bks, 26.Ot1, 27.Ot2,
%28.Edu, 29.Cha

ind = 34;

cons = [consdat(:,ind+1),...%foodhome
    consdat(:,ind+2)+consdat(:,ind+3),...%foodout + foodwork
    consdat(:,ind+4),...%tobacoo
    consdat(:,ind+5),...%alcohol
    consdat(:,ind+6),...%niteclub
    consdat(:,ind+7),...%clothes
    consdat(:,ind+8),...%tailor
    consdat(:,ind+9),...%jewelry
    consdat(:,ind+11),...%healthbeau
    consdat(:,ind+12),...%renthome
    consdat(:,ind+13),...%rentothr
    consdat(:,ind+14)+consdat(:,ind+15),...%furnish + housuppl
    consdat(:,ind+16)+consdat(:,ind+17)+consdat(:,ind+18)+consdat(:,ind+19),...%elect+gas+wager+homefuel
    consdat(:,ind+20),...%telephon
    consdat(:,ind+27),...%health ins
    consdat(:,ind+22)+consdat(:,ind+23)+consdat(:,ind+24)+consdat(:,ind+25)+consdat(:,ind+26),...%drugs+orthopd+doctors+hospital+nurshome
    consdat(:,ind+28),...%busisev
    consdat(:,ind+29),...%lifeins
    consdat(:,ind+30)+consdat(:,ind+31),...%autos+ parts
    consdat(:,ind+32),...%carservs
    consdat(:,ind+33),...%gasoline
    consdat(:,ind+35),...%autoins
    consdat(:,ind+36),...%masstran
    consdat(:,ind+37)+consdat(:,ind+38),...%othrtrans + airfare
    consdat(:,ind+39),...%books
    consdat(:,ind+40),...%pubs
    consdat(:,ind+41)+consdat(:,ind+42),...%
    consdat(:,ind+44)+consdat(:,ind+45)+consdat(:,ind+46),...%
    consdat(:,ind+47)];

%OLD VERSION BASED ON PDF DOCUMENTATION, CURRENT VERSION BASED ON LIBRARY LIST IN
%david_main.sas (library list downloaded from NBER website 
% cons = [consdat(:,ind+2),...
%     consdat(:,ind+1)+consdat(:,ind+3),...
%     consdat(:,ind+4),...
%     consdat(:,ind+6),...
%     consdat(:,ind+5),...
%     consdat(:,ind+7),...
%     consdat(:,ind+8),...
%     consdat(:,ind+9),...
%     consdat(:,ind+11),...
%     consdat(:,ind+12),...
%     consdat(:,ind+13),...
%     consdat(:,ind+14)+consdat(:,ind+15),...
%     consdat(:,ind+16)+consdat(:,ind+17)+consdat(:,ind+18)+consdat(:,ind+19),...
%     consdat(:,ind+20),...
%     consdat(:,ind+21),...
%     consdat(:,ind+22)+consdat(:,ind+23)+consdat(:,ind+24)+consdat(:,ind+25)+consdat(:,ind+26)+consdat(:,ind+27),...
%     consdat(:,ind+28),...
%     consdat(:,ind+29),...
%     consdat(:,ind+30),...
%     consdat(:,ind+31)+consdat(:,ind+32),...
%     consdat(:,ind+33),...
%     consdat(:,ind+35),...
%     consdat(:,ind+36)+consdat(:,ind+37),...
%     consdat(:,ind+38),...
%     consdat(:,ind+39)+consdat(:,ind+40),...
%     consdat(:,ind+41),...
%     consdat(:,ind+42),...
%     consdat(:,ind+44)+consdat(:,ind+45)+consdat(:,ind+46),...
%     consdat(:,ind+47)];

%family characteristics: [totwt,adjwt,over50,black,married,male,year]
%(WARNING: NEED TO CONFIRM THAT OVER50 IS ACTUALLY OVER40.  I THINK IT IS
%(I CHANGED THE COMPILATION PROGRAM)
char = [consdat(:,5),consdat(:,6),consdat(:,end-5:end-1)];

%sort mainspread
order = [9;6;3;8;19;2;26;1;28;30;5;4;12;27;25;20;10;11;23;24;21;14;29;7;13;16;17;22;15;18;31];
mainspread_sort = sortrows([order,mainspread],1);
mainspread_sort = mainspread_sort(:,2:end);

%price: same order as cons
price = mainspread_sort(1:29,5:end-2)';

%visibility indexes: column ordr the same, row order: Not black under 40,
%Not black over 40, Black under 40, Black over 40
vin = mainspread_sort(1:29,1:4)';
 
clearvars -except cons price vin char

cd '/gpfs/home/dcj138/work/ConspCons/'

save dat