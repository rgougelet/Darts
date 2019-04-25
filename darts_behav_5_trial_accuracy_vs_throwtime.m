% Initialize paths
clear; clc;
scriptdir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis';
cd(scriptdir)

% load XLSX SNAP data
[num,txt,raw] = xlsread('.\behavioral_data.xlsx');
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end
subj_ns = unique(xlsx.subject)';

% Trial level throwtime vs. distance
abs_dists = xlsx.distance(xlsx.distance < 6 & xlsx.throwtime < 5 & xlsx.pres == 0);
abs_throwtimes = xlsx.throwtime(xlsx.distance < 6 & xlsx.throwtime < 5 & xlsx.pres == 0);
pres_dists = xlsx.distance(xlsx.distance < 6 & xlsx.throwtime < 5 & xlsx.pres == 1);
pres_throwtimes = xlsx.throwtime(xlsx.distance < 6 & xlsx.throwtime < 5 & xlsx.pres == 1);

fig = figure('visible', 'off', 'color', 'w');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
mod = 1.5;
fig.PaperPosition = [0, 0, 8, 8]*mod;
fig.PaperSize = [8, 8]*mod;
fig.Position = [0.25, 0.25, 7.75, 7.75]*mod;
fig.Resize = 'off';
fig.InvertHardcopy = 'off';
title_mod = 1.3;

x= pres_throwtimes;
y = pres_dists;
pres_n = length(x);
[pres_r,p] = corr(x,y,'rows','complete')
s1 = scatter(x,y,10,'b','filled');

set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier', title_mod)
xlabel('Post-Cue Time Taken to Throw')
ylabel('Distance from Target')
l1 = lsline;

hold on
x= abs_throwtimes;
y = abs_dists;
abs_n = length(x);
[abs_r,p] = corr(x,y,'rows','complete')
s2 = scatter(x,y,10,'r','filled');
axis([0 6 0 6])
l2 = lsline;
l2(1).LineWidth = 4;
l2(1).Color = 'r';
l2(2).LineWidth = 4;
l2(2).Color = 'b';
legend([s1 s2],{'Planning','Memory'})
  
currentFigure = gcf;
title(currentFigure.Children(end), 'Trial-Level Accuracy vs. Throw Time');
print(fig,'-dsvg ','fig6.svg')