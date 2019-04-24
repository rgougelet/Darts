close all; clc; clear

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

mats = dir('*_aac.mat');

for mat_i = 1:length(mats)
	load(mats(mat_i).name)
	% Initialize paths
	scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
	cd(scriptdir)
	h = heatmap(aacs,'ColorbarVisible','off');
	h.XDisplayLabels = strings([1,length(h.XData)]);
	h.YDisplayLabels = strings([1,length(h.XData)]);
	
	set(gca,'fontsize', 20)
	% set(gca,'xtick',[])
	% ax = gca;
	% ax.XDisplayLabels = [];
	% ax.YDisplayLabels = [];
	%set(gca,'xticklabel',[])
	% set(gca,'TitleFontSizeMultiplier', title_mod)
	xlabel('Alpha Components')
	ylabel('Theta Components')
currentFigure = gcf;
title(currentFigure.Children(end), 'Theta to Alpha Amplitude-Amplitude Coupling');
print(fig,'-dsvg ',['./fig9_',num2str(mats(mat_i).name(1:end-4)),'.svg'])
end

