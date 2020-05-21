%% init
clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath(genpath([script_dir,'deps/']))
% data_dir = [script_dir,'data/Last Full Pipeline/'];
% data_dir = '/home/rgougelet/Desktop/darts/data/Redid the bandpass filters/';
addpath(genpath(data_dir))

subjs_to_include = {
	'571'
	'579'
	'580'
	'607'
	'608'
	'616'
	'619'
	'621'
	'627'
	'631'
	} ;
srate = 512 ;
load_r = load('lap_ic_r.mat') ;

% convert subject ids to numbers
for subj_i = 1:length(subjs_to_include)
	for row_i = 1:length(load_r.r)
		if load_r.r(row_i).subject == str2num(subjs_to_include{subj_i});
			load_r.r(row_i).subj_i = subj_i; 
		end
	end
end

% rename variables
r = struct2table(load_r.r);
r.throwtime = r.eeg_throwtime;
r.eeg_throwtime = [];
r.delaytime = r.eeg_delaytime;
r.eeg_delaytime = [];
r.delaystilltime = [];
load_r.r = table2struct(r);

% remove nans
r = struct2table(load_r.r);
r(any(isnan(r.Variables),2),:) = [] ; % remove all nans
r.distance = r.distance*6.55; % convert to real-life distance
load_r.r = table2struct(r);

%% PLOT - distance X subject X condition
r = struct2table(load_r.r) ;
clc; close all;
fig = figure('Position', [0 0 1366 786], 'Visible', 'on');
h =	tight_subplot(1,2,.045,.08);
t = annotation('textbox',[.5,.98,0,0],...
	'string','Dart''s Distance from Bull''s-Eye by Subject x Condition',...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);

% plot target absent
set(gcf,'CurrentAxes',h(1));
r = struct2table(load_r.r) ;
r(r.pres==1,:)=[]; % remove target present trials
for subj_i = 1:10
	v{subj_i} = r.distance(r.subj_i==subj_i);
end
violin(v, 'facecolor','b','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.distance, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',.75)
title('Target Absent'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Distance (cm)'); ylim([-1, 40]);
yline(mean(r.distance), 'LineWidth',2, 'Color','b')
legend({'Subject Distributions','Subject Means','Total Mean'})

% plot target present
set(gcf,'CurrentAxes',h(2));
r = struct2table(load_r.r) ;
r(r.pres==0,:)=[]; % remove target absent trials
for subj_i = 1:10
	v{subj_i} = r.distance(r.subj_i==subj_i);
end
violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.distance, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',.75)
title('Target Present'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Distance (cm)'); ylim([-1, 40]);
yline(mean(r.distance), 'LineWidth',2, 'Color','r')
legend({'Subject Distributions','Subject Means','Total Mean'})

tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% PLOT - throwtime X subject X condition
r = struct2table(load_r.r) ;
clc; close all;
fig = figure('Position', [0 0 1366 786], 'Visible', 'on');
h =	tight_subplot(1,2,.045,.1);
t = annotation('textbox',[.5,.98,0,0],...
	'string','Time Taken to Throw by Subject x Condition',...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);

% plot target absent
set(gcf,'CurrentAxes',h(1));
r = struct2table(load_r.r) ;
r(r.pres==1,:)=[]; % remove target present trials
for subj_i = 1:10
	v{subj_i} = r.throwtime(r.subj_i==subj_i);
end
violin(v, 'facecolor','b','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.throwtime, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',.75)
title('Target Absent'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Throw time (seconds)'); ylim([0.8, 3.5]);
yline(mean(r.throwtime), 'LineWidth',2, 'Color','b')
legend({'Subject Distributions','Subject Means','Total Mean'})

% plot target present
set(gcf,'CurrentAxes',h(2));
r = struct2table(load_r.r) ;
r(r.pres==0,:)=[]; % remove target absent trials
for subj_i = 1:10
	v{subj_i} = r.throwtime(r.subj_i==subj_i);
end
violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.throwtime, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',1)
yline(mean(r.throwtime), 'LineWidth',2, 'Color','r')
title('Target Present'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Throw time (seconds)'); ylim([0.8, 3.5]);
legend({'Subject Distributions','Subject Means','Total Mean'})

tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% PLOT - eeg vars X subject X condition
clc; close all;
for var = {'delay_Front_theta','delay_Front_alpha','delay_Front_gamma',...
	'pre_throw_Front_theta', 'pre_throw_Front_alpha', 'pre_throw_Front_gamma'}
	fig = figure('Position', [0 0 1366 786], 'Visible', 'on');
	h =	tight_subplot(1,2,.045,.1);
	t = annotation('textbox',[.5,.98,0,0],...
		'string',[(var{:}),' by Subject x Condition'],...
		'FitBoxToText','on',...
		'LineStyle','none',...
		'HorizontalAlignment','center',...
		'FontSize',14);
	
	% plot target absent
	set(gcf,'CurrentAxes',h(1));
	r = struct2table(load_r.r) ;
	r(r.pres==1,:)=[]; % remove target present trials
	for subj_i = 1:10
		v{subj_i} = r.(var{:})(r.subj_i==subj_i);
		% normalize each subject's distribution
% 		v{subj_i} = (v{subj_i}-mean(v{subj_i}))./std(v{subj_i});
	end
	violin(v, 'facecolor','b','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
	boxplot(r.(var{:}), r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
	set(findobj(gca,'type','line'),'linew',.75)
	title('Target Absent'); set(gca,'FontSize',12);
	xlabel('Subject'); ylabel('Throw time (seconds)');
	% ylim([0.8, 3.5]);
	yline(mean(r.(var{:})), 'LineWidth',2, 'Color','b')
	legend({'Subject Distributions','Subject Means','Total Mean'})
	
	% plot target present
	set(gcf,'CurrentAxes',h(2));
	r = struct2table(load_r.r) ;
	r(r.pres==0,:)=[]; % remove target absent trials
	for subj_i = 1:10
		v{subj_i} = r.(var{:})(r.subj_i==subj_i);
		% normalize each subject's distribution
% 		v{subj_i} = (v{subj_i}-mean(v{subj_i}))./std(v{subj_i});
	end
	violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
	boxplot(r.(var{:}), r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
	set(findobj(gca,'type','line'),'linew',1)
	yline(mean(r.(var{:})), 'LineWidth',2, 'Color','r')
	title('Target Present'); set(gca,'FontSize',12);
	xlabel('Subject'); ylabel('Throw time (seconds)');
	% ylim([0.8, 3.5]);
	legend({'Subject Distributions','Subject Means','Total Mean'})
	
	% tic; print(fig,t.String,'-dsvg','-r300'); toc;
end