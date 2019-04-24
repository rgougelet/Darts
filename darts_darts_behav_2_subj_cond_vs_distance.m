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

%% Generate within-subject stats
s_abs_mean_dists=[];s_pres_mean_dists=[];
s_abs_var_dists=[];s_pres_var_dists=[];
for subj_n = subj_ns
	subj = structfun(@(F) F([xlsx.subject]==subj_n), xlsx, 'uniform', 0);
	s_abs_mean_dists(end+1) = nanmean(subj.distance(subj.pres ==0));
	s_pres_mean_dists(end+1) = nanmean(subj.distance(subj.pres == 1));
	s_abs_var_dists(end+1) = nanvar(subj.distance(subj.pres == 0));
	s_pres_var_dists(end+1) = nanvar(subj.distance(subj.pres == 1));
end

%% Within-subject distance mean abs/pres
clc;
[H,P,CI,STATS] = ttest(s_abs_mean_dists-s_pres_mean_dists)

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

plot(s_abs_mean_dists,s_pres_mean_dists, 'o', 'MarkerSize',12, 'MarkerFaceColor','b'); hold on;
axis([1.3 2.2 1.3 2.2])
plot(0:0.01:3,0:0.01:3, 'Linewidth',3)
xlabel('Subject Average Distance from Target from Memory');
ylabel('Subject Average Distance from Target while Planning');
set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier', title_mod)

title_str = 'Within-Subject Differences in Accuracy';
title(title_str);
print(fig,'-dsvg ','fig2.svg')

%% Within-subject distance variance abs/pres

[H,P,CI,STATS] = ttest(s_abs_var_dists-s_pres_var_dists)
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

plot(s_abs_var_dists,s_pres_var_dists, 'o', 'MarkerSize',12, 'MarkerFaceColor','b'); hold on;
axis([0.5 2 0.5 2])
plot(0:0.01:3,0:0.01:3, 'Linewidth',3);
set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier', title_mod)
xlabel('Subject Variance Distance from Target from Memory');
ylabel('Subject Variance Distance from Target while Planning');

title_str = 'Within-Subject Differences in Precision';
title(title_str);
print(fig,'-dsvg ','fig3.svg')

%% trash
% c = categorical({'foo'});
% errorbar(mean(s_pres_mean_dists)-mean(s_abs_mean_dists),range(CI)/2,'.k','CapSize',100, 'Linewidth', 3)
% set(gca, 'XTickLabel', {''})
set(gca, 'YTick', [0 .1,.2,.3])
