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

% Trial-level throw times abs/ pres
abs_dists = xlsx.distance(xlsx.pres==0 & xlsx.distance < 6 & xlsx.throwtime < 5);
abs_delays = xlsx.delaystilltime(xlsx.pres==0 & xlsx.distance < 6 & xlsx.throwtime < 5);
abs_throwtimes = xlsx.throwtime(xlsx.pres==0 & xlsx.distance < 6 & xlsx.throwtime < 5);
abs_delay_tots = abs_delays+abs_throwtimes;

pres_dists = xlsx.distance(xlsx.pres==1& xlsx.distance < 6 & xlsx.throwtime < 5);
pres_delays = xlsx.delaystilltime(xlsx.pres==1& xlsx.distance < 6 & xlsx.throwtime < 5);
pres_throwtimes = xlsx.throwtime(xlsx.pres==1 & xlsx.distance < 6 & xlsx.throwtime < 5);
pres_delay_tots = pres_delays+pres_throwtimes;

% Histograms
close all
bins = 0:0.5:6;
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

sax1 = subplot(2,1,1);
histogram(abs_throwtimes,bins);
l11 = line([mean(abs_throwtimes), mean(abs_throwtimes)], [0 max(sax1.YLim)], 'color', 'r', 'Linewidth',3);
l21 = line([median(abs_throwtimes), median(abs_throwtimes)], [0 max(sax1.YLim)], 'color', 'g', 'Linewidth',3);
legend([l11,l21],'Mean','Median')
set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier', title_mod)
xlabel('Post-Cue Time Taken to Throw from \bfMemory')

sax2 = subplot(2,1,2);
histogram(pres_throwtimes,bins);
l12 = line([mean(pres_throwtimes), mean(pres_throwtimes)], [0 200], 'color', 'r', 'Linewidth',3);
l22 = line([median(pres_throwtimes), median(pres_throwtimes)], [0 200], 'color', 'g', 'Linewidth',3);
legend([l12,l22],'Mean','Median')
set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier', title_mod)
xlabel('Post-Cue Time Taken to Throw while \bfPlanning')

currentFigure = gcf;
title(currentFigure.Children(end), 'Time Taken to Throw');
print(fig,'-dsvg ','fig4.svg')

%% Stats
[H,P,CI,STATS] = ttest2(abs_throwtimes,pres_throwtimes)
[H,P,CI,STATS] = ttest2(abs_throwtimes,pres_throwtimes, 'vartype', 'unequal')
[H,P,STATS] = ranksum(abs_throwtimes,pres_throwtimes)
dists = [abs_throwtimes,padarray(pres_throwtimes,length(abs_throwtimes)-length(pres_throwtimes), NaN,'post')];
% p = vartestn(dists,'TestType','LeveneAbsolute')