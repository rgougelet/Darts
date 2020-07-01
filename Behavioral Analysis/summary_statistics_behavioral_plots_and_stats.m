
clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);

addpath([script_dir,'eeglab/']);
addpath(genpath([script_dir,'deps/']));
data_dir = [script_dir,'data/'];
addpath(genpath(data_dir));

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

% load XLSX SNAP data
[num,txt,~] = xlsread([data_dir,'behavioral_data_reduced.xlsx']);
r = struct;
headers = txt(1,:);
for k=1:numel(headers)
	for ri = 1:length(num)
		r(ri).(headers{k})=num(ri,k);
	end
end

load_r.r = r;
% convert subject ids to numbers
for subj_i = 1:length(subjs_to_include)
	for row_i = 1:length(load_r.r)
		if load_r.r(row_i).subject == str2num(subjs_to_include{subj_i}); load_r.r(row_i).subject = subj_i; end
	end
end
% rename variables
% r = struct2table(load_r.r);
% load_r.r = table2struct(r);

% remove nans
r = struct2table(load_r.r);
r(any(isnan(r.Variables),2),:) = [] ; % remove all nans
r.distance = r.distance*6.55; % convert to real-life distance
load_r.r = table2struct(r);

fig_dir = [script_dir,'figures/'];
format short g

%% STAT - summary statistics distance X condition - sig
r = struct2table(load_r.r) ;
abs_subj_means_distance = []; pres_subj_means_distance = [];
for subj_i = 1:length(subjs_to_include)
	abs_subj_means_distance(end+1,:) = mean(r.distance(r.subject == subj_i & r.pres == 0));
	pres_subj_means_distance(end+1,:) = mean(r.distance(r.subject == subj_i & r.pres == 1));
end
% [~,p,~,stats] = ttest(abs_subj_means_distance,pres_subj_means_distance)
[~,p,~,se] = perm_ttest(abs_subj_means_distance,pres_subj_means_distance, 100000)
% [~,p] = ttest2(abs_subj_means_distance,pres_subj_means_distance)
[~,p] = perm_ttest2(abs_subj_means_distance,pres_subj_means_distance,'nperm', 10000)
%% PLOT - summary statistics distances X subject X condition
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
	v{subj_i} = r.distance(r.subject==subj_i);
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
	v{subj_i} = r.distance(r.subject==subj_i);
end
violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.distance, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',.75)
title('Target Present'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Distance (cm)'); ylim([-1, 40]);
yline(mean(r.distance), 'LineWidth',2, 'Color','r')
legend({'Subject Distributions','Subject Means','Total Mean'})

tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% STAT - summary statistics throwtime X condition - sig
r = struct2table(load_r.r) ;
abs_subj_means_throwtime = []; pres_subj_means_throwtime = [];
for subj_i = 1:length(subjs_to_include)
	abs_subj_means_throwtime(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.pres == 0));
	pres_subj_means_throwtime(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.pres == 1));
end
% [h,p,~,stats] = ttest(abs_subj_means_distance,pres_subj_means_distance)
[xbar,p,~,se] = perm_ttest(abs_subj_means_throwtime,pres_subj_means_throwtime, 10000);
% stats.sd/sqrt(length(abs_subj_means_distance-pres_subj_means_distance))
% [h,p] = ttest2(abs_subj_means_distance,pres_subj_means_distance)
% [h,p] = perm_ttest2(abs_subj_means_distance,pres_subj_means_distance,'nperm', 10000)
%% PLOT - summary statistics throwtime X subject X condition
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
	v{subj_i} = r.throwtime(r.subject==subj_i);
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
	v{subj_i} = r.throwtime(r.subject==subj_i);
end
violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.throwtime, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
set(findobj(gca,'type','line'),'linew',1)
yline(mean(r.throwtime), 'LineWidth',2, 'Color','r')
title('Target Present'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Throw time (seconds)'); ylim([0.8, 3.5]);
legend({'Subject Distributions','Subject Means','Total Mean'})

tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% STAT - summary statistic subject means for each delay correlation against distance - near sig
r = struct2table(load_r.r);
abs_subj_means_per_delay = []; pres_subj_means_per_delay = [];
for subj_i = 1:length(subjs_to_include)
	for delay = 3:9
		abs_subj_means_per_delay(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		pres_subj_means_per_delay(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==1));
	end
end

ds = repmat(3:9,10,1);
ds = ds(:);
as = abs_subj_means_per_delay(:);
ps = pres_subj_means_per_delay(:);

% [rho, p] = corr(ds,as); [rho, p]
% [rho, p] = corr(ds,as, 'type','Spearman'); [rho, p]
% [rho, p] = corr(ds,as, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(ds,as); [rho, p]
% [rho, p] = studentized_perm_corr(ds,as); [rho, p]

% [rho, p] = corr(ds,ps); [rho, p]
% [rho, p] = corr(ds,ps, 'type','Spearman'); [rho, p]
% [rho, p] = corr(ds,ps, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(ds,ps); [rho, p]
% [rho, p] = studentized_perm_corr(ds,ps); [rho, p]

r1 = corr(ds,as);
n1 = length(as);
r2 = corr(ds,ps);
n2 = length(ps);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

% [rho, p] = corr(as-ps,ds); [rho, p]
% [rho, p] = corr(as-ps,ds, 'type','Spearman'); [rho, p]
% [rho, p] = corr(as-ps,ds, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(as-ps,ds); [rho, p]
% [rho, p] = studentized_perm_corr(as-ps,ds); [rho, p]
%% PLOT - summary statistic subject means for each delay correlation against distance
r = struct2table(load_r.r);
abs_subj_means_per_delay = []; pres_subj_means_per_delay = [];
for subj_i = 1:length(subjs_to_include)
	for delay = 3:9
		abs_subj_means_per_delay(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		pres_subj_means_per_delay(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==1));
	end
end

ds = repmat(3:9,10,1);
ds = ds(:);
as = abs_subj_means_per_delay(:);
ps = pres_subj_means_per_delay(:);
figure; plot(ds,as,'o')
lsline
figure; plot(3:9, abs_subj_means_per_delay,'o')
lsline
plot(3:9, pres_subj_means_per_delay,'o')

%% STAT - summary statistic subject means for throwtime against distance - near sig
r = struct2table(load_r.r);
abs_subj_means_per_delay = []; pres_subj_means_per_delay = [];
for subj_i = 1:length(subjs_to_include)
	for delay = 3:9
		abs_subj_means_per_delay(subj_i,delay-2) = mean(r.throwtime(r.subject==subj_i & r.delay==delay & r.pres==0));
		pres_subj_means_per_delay(subj_i,delay-2) = mean(r.throwtime(r.subject==subj_i & r.delay==delay & r.pres==1));
	end
end

ds = repmat(3:9,10,1);
ds = ds(:);
as = abs_subj_means_per_delay(:);
ps = pres_subj_means_per_delay(:);

% [rho, p] = corr(ds,as); [rho, p]
% [rho, p] = corr(ds,as, 'type','Spearman'); [rho, p]
% [rho, p] = corr(ds,as, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(ds,as); [rho, p]
% [rho, p] = studentized_perm_corr(ds,as); [rho, p]

% [rho, p] = corr(ds,ps); [rho, p]
% [rho, p] = corr(ds,ps, 'type','Spearman'); [rho, p]
% [rho, p] = corr(ds,ps, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(ds,ps); [rho, p]
% [rho, p] = studentized_perm_corr(ds,ps); [rho, p]

r1 = corr(ds,as);
n1 = length(as);
r2 = corr(ds,ps);
n2 = length(ps);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

% [rho, p] = corr(as-ps,ds); [rho, p]
% [rho, p] = corr(as-ps,ds, 'type','Spearman'); [rho, p]
% [rho, p] = corr(as-ps,ds, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(as-ps,ds, 100000); [rho, p]
[rho, p] = studentized_perm_corr(as-ps,ds); [rho, p]


%% dont bother with below?
%% STAT - summary statistic delay vs mean distance per delay - near sig, not really
r = struct2table(load_r.r) ;

abs_delay_distance_means = [];
pres_delay_distance_means = [];
% average the subject means within each second of delay
for delay = 3:9
	abs_delay_distance_means(end+1,:) = mean(r.distance(r.delay==delay& r.pres==0));
	pres_delay_distance_means(end+1,:) = mean(r.distance(r.delay==delay& r.pres==1));
end
disp('Summary Statistics')
disp('Target Absent')
[rho, p] = corr((3:9)',abs_delay_distance_means); [rho, p]
[rho, p] = corr((3:9)',abs_delay_distance_means, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',abs_delay_distance_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',abs_delay_distance_means); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',abs_delay_distance_means); [rho, p]
disp('Target Present')
[rho, p] = corr((3:9)',pres_delay_distance_means); [rho, p]
[rho, p] = corr((3:9)',pres_delay_distance_means, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',pres_delay_distance_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',pres_delay_distance_means); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',pres_delay_distance_means); [rho, p]
disp('Target Absent - Target Present')
d = abs_delay_distance_means-pres_delay_distance_means;
[rho, p] = perm_corr((3:9)',d); [rho, p]
[rho, p] = corr((3:9)',d); [rho, p]
[rho, p] = corr((3:9)',d, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',d, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',d); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',d); [rho, p]

r1 = corr((3:9)',abs_delay_distance_means);
n1 = length(abs_delay_distance_means);
r2 = corr((3:9)',pres_delay_distance_means);
n2 = length(pres_delay_distance_means);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)
%% PLOT - summary statistic delay vs mean distance per delay
close all; figure;
subplot(2,1,1)
plot(3:9,abs_delay_distance_means,'ko')
xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Delay (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

subplot(2,1,2)
plot(3:9,pres_delay_distance_means,'ko')
xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Delay (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)


%% STAT - summary statistic distance vs throwtime correlation using subject means - near sig 
r = struct2table(load_r.r) ;
abs_subj_means_distance = []; abs_subj_means_throwtime = [];
pres_subj_means_distance = []; pres_subj_means_throwtime = [];
for subj_i = 1:length(subjs_to_include)
		abs_subj_means_throwtime(end+1,:) = mean(r.throwtime(r.subject==subj_i & r.pres==0));
		abs_subj_means_distance(end+1,:) = mean(r.distance(r.subject==subj_i & r.pres==0));
		pres_subj_means_throwtime(end+1,:) = mean(r.throwtime(r.subject==subj_i & r.pres==1));
		pres_subj_means_distance(end+1,:) = mean(r.distance(r.subject==subj_i & r.pres==1));
end
disp('Summary Statistics')
disp('Target Absent')
[rho, p] = corr(abs_subj_means_throwtime,abs_subj_means_distance); [rho, p]
[rho, p] = corr(abs_subj_means_throwtime,abs_subj_means_distance, 'type','Spearman'); [rho, p]
[rho, p] = corr(abs_subj_means_throwtime,abs_subj_means_distance, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(abs_subj_means_throwtime,abs_subj_means_distance); [rho, p]
[rho, p] = studentized_perm_corr(abs_subj_means_throwtime,abs_subj_means_distance); [rho, p]

disp('Target Present')
[rho, p] = corr(pres_subj_means_throwtime,pres_subj_means_distance); [rho, p]
[rho, p] = corr(pres_subj_means_throwtime,pres_subj_means_distance, 'type','Spearman'); [rho, p]
[rho, p] = corr(pres_subj_means_throwtime,pres_subj_means_distance, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(pres_subj_means_throwtime,pres_subj_means_distance); [rho, p]
[rho, p] = studentized_perm_corr(pres_subj_means_throwtime,pres_subj_means_distance); [rho, p]

r1 = corr(abs_subj_means_throwtime,abs_subj_means_distance);
n1 = length(abs_subj_means_throwtime);
r2 = corr(pres_subj_means_throwtime,pres_subj_means_distance);
n2 = length(pres_subj_means_throwtime);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

[rho, p] = corr(abs_subj_means_throwtime-pres_subj_means_throwtime,abs_subj_means_distance-pres_subj_means_distance); [rho, p]
[rho, p] = corr(abs_subj_means_throwtime-pres_subj_means_throwtime,abs_subj_means_distance-pres_subj_means_distance, 'type','Spearman'); [rho, p]
[rho, p] = corr(abs_subj_means_throwtime-pres_subj_means_throwtime,abs_subj_means_distance-pres_subj_means_distance, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(abs_subj_means_throwtime-pres_subj_means_throwtime,abs_subj_means_distance-pres_subj_means_distance); [rho, p]
[rho, p] = studentized_perm_corr(abs_subj_means_throwtime-pres_subj_means_throwtime,abs_subj_means_distance-pres_subj_means_distance); [rho, p]
%% PLOT - summary statistic distance vs throwtime correlation using subject means
close all;
fig = figure('Position', [0 0 1366 786], 'Visible', 'on'); 
h =	tight_subplot(1,2,.045,.1);
t = annotation('textbox',[.5,.98,0,0],...
	'string','Subject Average Time Taken to Throw vs Subject Average Distance from Target',...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);
subplot(3,1,1)
plot(abs_subj_means_throwtime,abs_subj_means_distance,'r*')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
title('Target Absent');
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

subplot(3,1,2)
plot(pres_subj_means_throwtime,pres_subj_means_distance,'r*')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
title('Target Present');
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

subplot(3,1,3)
plot(abs_subj_means_throwtime-pres_subj_means_throwtime,abs_subj_means_distance-pres_subj_means_distance,'r*')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
title('Target Absent - Target Present');
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

% %% STAT - summary statistic distance vs throwtime correlation using subject x delay means
% r = struct2table(load_r.r) ;
% abs_cond_subj_throwtimes= []; abs_subj_means_throwtimes_conds = [];
% pres_cond_subj_throwtimes = []; pres_subj_means_throwtimes_conds = [];
% for subj_i = 1:10
% 	for delay = 3:9
% 		abs_cond_subj_throwtimes(subj_i,delay) = mean(r.throwtime(r.subject == subj_i & r.delay==delay));
% 		pres_cond_subj_throwtimes(subj_i,delay) = mean(r.throwtime(r.subject == subj_i & r.delay==delay));
% 	end
% end
% disp('Summary Statistics')
% disp('Target Absent')
% [rho, p] = corr((3:9)',abs_subj_means_throwtimes_conds); [rho, p]
% [rho, p] = corr((3:9)',abs_subj_means_throwtimes_conds, 'type','Spearman'); [rho, p]
% [rho, p] = corr((3:9)',abs_subj_means_throwtimes_conds, 'type','Kendall'); [rho, p]
% [rho,p] = perm_corr((3:9)',abs_subj_means_throwtimes_conds); [rho, p]
% [rho,p] = studentized_perm_corr((3:9)',abs_subj_means_throwtimes_conds); [rho, p]
% disp('Target Present')
% [rho, p] = corr((3:9)',pres_subj_means_throwtimes_conds); [rho, p]
% [rho, p] = corr((3:9)',pres_subj_means_throwtimes_conds, 'type','Spearman'); [rho, p]
% [rho, p] = corr((3:9)',pres_subj_means_throwtimes_conds, 'type','Kendall'); [rho, p]
% [rho,p] = perm_corr((3:9)',pres_subj_means_throwtimes_conds); [rho, p]
% [rho,p] = studentized_perm_corr((3:9)',pres_subj_means_throwtimes_conds); [rho, p]
% disp('Target Absent - Target Present')
% d = abs_subj_means_throwtimes_conds-pres_subj_means_throwtimes_conds;
% [rho, p] = perm_corr((3:9)',d); [rho, p]
% [rho, p] = corr((3:9)',d); [rho, p]
% [rho, p] = corr((3:9)',d, 'type','Spearman'); [rho, p]
% [rho, p] = corr((3:9)',d, 'type','Kendall'); [rho, p]
% [rho,p] = perm_corr((3:9)',d); [rho, p]
% [rho,p] = studentized_perm_corr((3:9)',d); [rho, p]
% 
% r1 = corr((3:9)',abs_subj_means_distance_conds);
% n1 = length(abs_subj_means_distance_conds);
% r2 = corr((3:9)',pres_subj_means_distance_conds);
% n2 = length(pres_subj_means_distance_conds);
% [p,z] = corr_rtest(r1,r2,n1,n2); p(2)
% %% PLOT -  summary statistic distance vs throwtime correlation using subject x condition means
% subplot(2,1,1)
% plot(3:9,abs_subj_means_distance_conds,'ko')
% % xlim([2.5, 9.5])
% % ylim([0.75, 7])
% xlabel('Throw time (seconds)')
% ylabel('Distance (cm)')
% lsline;
% set(gca,'FontSize',12)
% 
% subplot(2,1,2)
% plot(3:9,pres_subj_means_distance_conds,'ko')
% % xlim([2.5, 9.5])
% % ylim([0.75, 7])
% xlabel('Throw time (seconds)')
% ylabel('Distance (cm)')
% lsline;
% set(gca,'FontSize',12)






















