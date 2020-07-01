%% init
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

%% todo: verify/remove sqrts
% dv distance = idv delay, idv condition, dv throwtime, i/dv theta, i/dv alpha, (1+subject|throwtime)
% dv throwtime = idv delay, idv condition, i/dv theta, i/dv alpha, (1+subject)
% dv distance = 1 x 2 x 1 x 1 x 1
% dv throwtime = 1 x 2 x 1 x 1

% done, distance by subject per condition
% check theta variability, distance vs theta per condition
% check alpha variability, distance vs alpha per condition

%% PLOT - single-trial distance X condition
clc; close all;
r = struct2table(load_r.r) ;
fig = figure('Position', [0 0 1366 786], 'Visible', 'on'); 
h1 = histogram(r.distance(~r.pres,:), 'BinEdges',0:1:49, 'FaceAlpha',0.15); hold on;
h2 = histogram(r.distance(~~r.pres,:), 'BinEdges',0:1:49, 'FaceAlpha',0.15);hold on;
o = get(h1.Parent, 'ColorOrder');
h1v = h1.Values;
h2v = h2.Values; clf;
b = bar(.5:1:49,[h2v;h1v]', 1,'stacked', 'FaceAlpha',0.15);
b(1).FaceColor = o(2,:);
b(2).FaceColor = o(1,:);
ylabel('Trials');
xlabel('Distance (cm)');
xline(mean(r.distance(~~r.pres,:)), 'LineWidth',2, 'Color','r');
xline(mean(r.distance(~r.pres,:)), 'LineWidth',2, 'Color','b');
set(gca,'FontSize',12);
t = title('Dart''s Distance from Bull''s-Eye by Condition');
l = legend({'Target Present','Target Absent', 'Target Present Mean', 'Target Absent Mean'});
set(l,'position',[0.75 0.8 0.1 0.1]);
tic; print(fig,[fig_dir,t.String],'-dsvg','-r300'); toc;
%% STAT - sig - single-trial distance X condition
r = struct2table(load_r.r) ;
ta = r.distance(r.pres == 0);
tp = r.distance(r.pres == 1);
[h,p] = ttest2(ta,tp)
[h,p] = perm_ttest2(ta,tp,'nperm', 10000)

%% PLOT - single-trial throwtime X condition
clc; close all;
r = struct2table(load_r.r) ;
fig = figure('Position', [0 0 1366 786], 'Visible','on');
h1 = histogram(r.throwtime(~r.pres,:), 'BinEdges',0.8:.2:3.6, 'FaceAlpha',0.15); hold on;
h2 = histogram(r.throwtime(~~r.pres,:), 'BinEdges',0.8:.2:3.6, 'FaceAlpha',0.15); hold on;
o = get(h1.Parent, 'ColorOrder');
h1v = h1.Values;
h2v = h2.Values; clf;
b = bar(.9:.2:3.6,[h2v;h1v]', 1,'stacked', 'FaceAlpha',0.15);
b(1).FaceColor = o(2,:);
b(2).FaceColor = o(1,:);
ylabel('Trials');
xlabel('Throw time (seconds)');
xline(mean(r.throwtime(~~r.pres,:)), 'LineWidth',2, 'Color','r');
xline(mean(r.throwtime(~r.pres,:)), 'LineWidth',2, 'Color','b');
set(gca,'FontSize',12);
t = title('Time Taken to Throw by Condition');
l = legend({'Target Present','Target Absent', 'Target Present Mean', 'Target Absent Mean'});
set(l,'position',[0.75 0.8 0.1 0.1])
xlim([0.75, 3.5]);
tic; print(fig,[fig_dir,t.String],'-dsvg','-r300'); toc;
%% STAT - sig - single-trial throwtime X condition
r = struct2table(load_r.r) ;
ta = r.throwtime(r.pres == 0);
tp = r.throwtime(r.pres == 1);
[h,p] = ttest2(ta,tp)
[h,p] = perm_ttest2(ta,tp,'nperm', 100000)

%% PLOT - single-trial distance vs delay correlation by condition
close all;
fig = figure('Position', [0 0 1366 786], 'Visible','on');
t = annotation('textbox',[.5,.98,0,0],...
	'string','Dart''s Distance from Bull''s-Eye by Delay Time',...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);
r = struct2table(load_r.r) ;
abs_delay = r.delay(r.pres==0);
pres_delay = r.delay(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);

% single trial plot
subplot(2,1,1)
plot(abs_delay, sqrt(abs_distance),'b.'); hold on;
title('Target Absent');
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay Time (seconds)')
ylabel('Distance (sqrt(cm))')
l_abs = lsline();
l_abs.Color = 'b';
set(gca,'FontSize',12)

subplot(2,1,2)
plot(pres_delay, sqrt(pres_distance),'r.'); hold on;
title('Target Present');
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay Time (seconds)')
ylabel('Distance (sqrt(cm))')
l_pres = lsline();
l_pres.Color = 'r';
set(gca,'FontSize',12)
tic; print(fig,[fig_dir,t.String],'-dsvg','-r300'); toc;
%% STAT - near sig - single-trial distance vs delay correlation by condition
r = struct2table(load_r.r) ;
% r(1:10) = [];
abs_delay = r.delay(r.pres==0);
pres_delay = r.delay(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);
disp('Single Trial')
disp('Target Absent')
% [rho, p] = corr(abs_delay,abs_distance); [rho, p]
% [rho, p] = corr(abs_delay,abs_distance, 'type','Spearman'); [rho, p]
% [rho, p] = corr(abs_delay,abs_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(abs_delay,abs_distance); [rho, p]
% [rho, p,~,~] = studentized_perm_corr(abs_delay,abs_distance); [rho, p]
disp('Target Present')
% [rho, p] = corr(pres_delay,pres_distance); [rho, p]
% [rho, p] = corr(pres_delay,pres_distance, 'type','Spearman'); [rho, p]
% [rho, p] = corr(pres_delay,pres_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(pres_delay,pres_distance); [rho, p]
% [rho, p,~,~] = studentized_perm_corr(pres_delay,pres_distance); [rho, p]
r1 = corr(abs_delay,abs_distance);
n1 = length(abs_delay);
r2 = corr(pres_delay,pres_distance);
n2 = length(pres_delay);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

[rho, p,~,~] = perm_corr(r.delay,r.distance)

%% PLOT - single-trial distance vs throw time correlation by condition
clc; close all;
fig = figure('Position', [0 0 1366 786], 'Visible','on');
t = annotation('textbox',[.5,.98,0,0],...
	'string','Dart''s Distance from Bull''s-Eye by Throw Time',...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);
r = struct2table(load_r.r) ;
r(r.throwtime > 4,:) = [];
abs_throwtime = r.throwtime(r.pres==0);
pres_throwtime = r.throwtime(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);

% single trial plot
subplot(2,1,1)
plot(abs_throwtime, abs_distance,'b.'); hold on;
title('Target Absent');
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
l_abs = lsline();
l_abs.Color = 'b';
set(gca,'FontSize',12)

subplot(2,1,2)
plot(pres_throwtime, pres_distance,'r.'); hold on;
title('Target Present');
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
l_pres = lsline();
l_pres.Color = 'r';
set(gca,'FontSize',12)
tic; print(fig,t.String,'-dsvg','-r300'); toc;
%% STAT - sig - statistics single-trial distance vs throwtime correlation by condition
r = struct2table(load_r.r) ;
r(r.throwtime > 4,:) = [];
abs_throwtime = r.throwtime(r.pres==0);
pres_throwtime = r.throwtime(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);
% disp('Single Trial')
% disp('Target Absent')
% [rho, p] = corr(abs_throwtime,abs_distance); [rho, p]
% [rho, p] = corr(abs_throwtime,abs_distance, 'type','Spearman'); [rho, p]
% [rho, p] = corr(abs_throwtime,abs_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(abs_throwtime,abs_distance); [rho, p] % .012
% [rho, p,~,~] = studentized_perm_corr(abs_throwtime,abs_distance); [rho, p]
% disp('Target Present')
% [rho, p] = corr(pres_throwtime,pres_distance); [rho, p]
% [rho, p] = corr(pres_throwtime,pres_distance, 'type','Spearman'); [rho, p]
% [rho, p] = corr(pres_throwtime,pres_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(pres_throwtime,pres_distance); [rho, p]
% [rho, p,~,~] = studentized_perm_corr(pres_throwtime,pres_distance); [rho, p]
% r1 = corr(abs_throwtime,abs_distance);
% n1 = length(abs_throwtime);
% r2 = corr(pres_throwtime,pres_distance);
% n2 = length(pres_throwtime);
% [p,z] = corr_rtest(r1,r2,n1,n2); p(2)

%% single-trial regression
r = struct2table(load_r.r) ;
r.pres = categorical(r.pres);

% lme = fitlm(r,'distance ~ throwtime*delay*pres')
lme = fitlme(r,'distance ~ throwtime*delay*pres + (1+delay|subject)+ (1+throwtime|subject)', 'FitMethod', 'REML', 'StartMethod','random')
lme = fitlme(r,'distance ~ throwtime*delay*pres + (1+throwtime|subject)', 'FitMethod', 'REML', 'StartMethod','random')










