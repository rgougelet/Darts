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

% Subject-Level Accuracy vs Throw Time
s_abs_throwtimes = []; s_abs_dists = []; s_pres_throwtimes = []; s_pres_dists = [];
for subj_n = subj_ns
	subj = structfun(@(F) F([xlsx.subject]==subj_n), xlsx, 'uniform', 0);
	s_abs_throwtimes(end+1) = nanmean(subj.throwtime(subj.pres == 0& subj.distance < 6 & subj.throwtime < 6));
	s_abs_dists(end+1) = nanmean(subj.distance(subj.pres == 0& subj.distance < 6 & subj.throwtime < 6));
	s_pres_throwtimes(end+1) = nanmean(subj.throwtime(subj.pres == 1& subj.distance < 6 & subj.throwtime < 6));
	s_pres_dists(end+1) = nanmean(subj.distance(subj.pres == 1& subj.distance < 6 & subj.throwtime < 6));
end

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

x= s_pres_throwtimes';
y = s_pres_dists';
[pres_r,p] = corr(x,y,'rows','complete')
s1 = scatter(x,y,120,'b','filled');
axis([1 2.5 1 2.5])
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
dx = 0.01; dy = 0.01;
text(x+dx, y+dy, c,'b','fontsize', 20);
set(gca,'fontsize', 20)
set(gca,'TitleFontSizeMultiplier', title_mod)
xlabel('Post-Cue Time Taken to Throw ')
ylabel('Distance from Target')
l1 = lsline;

hold on
x= s_abs_throwtimes';
y = s_abs_dists';
s2 = scatter(x,y,120,'r','filled');
[abs_r,p] = corr(x,y,'rows','complete')
l2 = lsline;
l2(1).LineWidth = 4;
l2(1).Color = 'r';
l2(2).LineWidth = 4;
l2(2).Color = 'b';
a = [1:length(x)]'; b = num2str(a); c = cellstr(b);
dx = 0.01; dy = 0.01;
text(x+dx, y+dy, c,'r','fontsize', 20);

legend([s1 s2],{'Planning', 'Memory'})

currentFigure = gcf;
title(currentFigure.Children(end), 'Subject-Level Accuracy vs. Throw Time');
print(fig,'-dsvg ','fig5.svg')