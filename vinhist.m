clear;
load 4-7-12.mat;

ot = histc(ot_ind(type==2),1:29);
ot = ot/sum(ot);

vinmatch = figure(1);

lab = {'FdH',' FdO', 'Cig', 'AlH', 'AlO', 'Clo', 'Lry', 'Jwl', 'Brb','Hom', 'Htl', 'Fur', 'Utl', 'Tel', 'HIn', 'Med', 'Fee', 'LIn', 'Car', 'CMn', 'Gas', 'CIn', 'Bus', 'Air', 'Bks', 'Ot1', 'Ot2', 'Edu', 'Cha'};
lscatter(log(v{2}),log(ot),lab,'TextColor','black','FontSize',8);
xlim([-7 -1]);
ylim([-7 -1]);
hold on
plot(-7:-1,-7:-1,':r');
hold off
xlabel('Vindex Probabilities');
ylabel('Estimated Observation Type Probabilities');
print(vinmatch,'-dpdf','figures/vinmatch');

