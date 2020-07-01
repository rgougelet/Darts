function plot_abs_vs_pres_subjs(load_r, varname, tit, ylab,visible)
clc; close all;
fig = figure('Position', [0 0 1366 786], 'Visible', visible); 
h =	tight_subplot(1,2,.045,.08);
t = annotation('textbox',[.5,.98,0,0],...
	'string',[tit, ' by Subject x Condition'],...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);

% plot target absent
set(gcf,'CurrentAxes',h(1));
r = load_r;
r(r.TargetPresent,:)=[]; % remove target present trials
for subj_i = 1:10
	v{subj_i} = r.(varname)(r.Subject==subj_i);
end
violin(v, 'facecolor','b','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.(varname), r.SubjectID, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',.75)
title('Target Absent'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel(ylab); 
yline(mean(r.(varname)), 'LineWidth',2, 'Color','b');
legend({'Subject Distributions','Subject Means','Total Mean'}, 'Location','northwest')

% plot target present
set(gcf,'CurrentAxes',h(2));
r = load_r;
r(r.TargetAbsent,:)=[]; % remove target absent trials
for subj_i = 1:10
	v{subj_i} = r.(varname)(r.Subject==subj_i);
end
violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.(varname), r.SubjectID, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',.75)
title('Target Present'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel(ylab); 
yline(mean(r.(varname)), 'LineWidth',2, 'Color','r');
legend({'Subject Distributions','Subject Means','Total Mean'}, 'Location','northwest')
tic; print(fig,t.String,'-dsvg','-r300'); toc;