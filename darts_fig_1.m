% Initialize paths
clear;
scriptdir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis';
cd(scriptdir)

% load XLSX SNAP data
[num,txt,raw] = xlsread('.\behavioral_data.xlsx');
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end
subj_ns = unique(xlsx.subject)';

%% Trial level cond vs. distance
abs_dists = xlsx.distance(xlsx.pres==0 & xlsx.distance < 6);
pres_dists = xlsx.distance(xlsx.pres==1& xlsx.distance < 6);

% Histograms
close all
bins = 0:0.5:6;
fig1 = figure('visible', 'off', 'color', 'w');
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
mod = 1.5;
fig1.PaperPosition = [0, 0, 8, 8]*mod;
fig1.PaperSize = [8, 8]*mod;
fig1.Position = [0.25, 0.25, 7.75, 7.75]*mod;
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';
title_mod = 1.3;

sax1 = subplot(2,1,1);
histogram(abs_dists,bins);
l1 = line([mean(abs_dists), mean(abs_dists)], [0 max(sax1.YLim)], 'color', 'r', 'Linewidth',3);
l2 = line([median(abs_dists), median(abs_dists)], [0 max(sax1.YLim)], 'color', 'g', 'Linewidth',3);
legend([l1,l2],'Mean','Median')
set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier',title_mod)
xlabel('Distance (cm) from Target from \bfMemory')

sax2 = subplot(2,1,2);
histogram(pres_dists,bins);
l1 = line([mean(pres_dists), mean(pres_dists)], [0 max(sax2.YLim)], 'color', 'r', 'Linewidth',3);
l2 = line([median(pres_dists), median(pres_dists)], [0 max(sax2.YLim)], 'color', 'g', 'Linewidth',3);
legend([l1,l2],'Mean','Median')
set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier', title_mod)
xlabel('Distance (cm) from Target with \bfPlanning')

currentFigure = gcf;
title(currentFigure.Children(end), 'Trial Accuracy');
print(fig1,'-dsvg ',[scriptdir,'/fig1.svg'])

%% Stats
[H,P,CI,STATS] = ttest2(abs_dists,pres_dists)
[H,P,CI,STATS] = ttest2(abs_dists,pres_dists, 'vartype', 'unequal')
dists = [abs_dists,padarray(pres_dists,length(abs_dists)-length(pres_dists), NaN,'post')];
p = vartestn(dists,'TestType','Bartlett')