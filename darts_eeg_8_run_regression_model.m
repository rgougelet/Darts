%% init
clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
addpath([script_dir,'eeglab/']);
addpath([script_dir,'deps/']);
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
r = load_r.r;
mperm_corr = memoize(@perm_corr);

% convert subject ids to numbers
for subj_i = 1:length(subjs_to_include)
	for row_i = 1:length(load_r.r)
		if load_r.r(row_i).subject == str2num(subjs_to_include{subj_i}); load_r.r(row_i).subject = subj_i; end
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

%% plot distance X subject X condition
clc; close all;
fig = figure('Position', [0 0 1366 786], 'Visible', 'off'); 
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
title('Target Absent'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Distance (cm)'); ylim([-1, 40]);
set(findobj(gca,'type','line'),'linew',.75)

% plot target present
set(gcf,'CurrentAxes',h(2));
r = struct2table(load_r.r) ;
r(r.pres==0,:)=[]; % remove target absent trials
for subj_i = 1:10
	v{subj_i} = r.distance(r.subject==subj_i);
end
violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.distance, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
title('Target Present'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Distance (cm)'); ylim([-1, 40]);
set(findobj(gca,'type','line'),'linew',.75)

tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% plot throwtime X subject X condition
clc; close all;
fig = figure('Position', [0 0 1366 786], 'Visible', 'off'); 
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
title('Target Absent'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Throw time (seconds)'); ylim([0.8, 3.5]);
set(findobj(gca,'type','line'),'linew',.75)

% plot target present
set(gcf,'CurrentAxes',h(2));
r = struct2table(load_r.r) ;
r(r.pres==0,:)=[]; % remove target absent trials
for subj_i = 1:10
	v{subj_i} = r.throwtime(r.subject==subj_i);
end
violin(v, 'facecolor','r','facealpha',0.25, 'mc','k','medc','','edgecolor',[]);
boxplot(r.throwtime, r.subject, 'widths',0.15, 'colors','k', 'Symbol','+','OutlierSize',3); hold on;
title('Target Present'); set(gca,'FontSize',12);
xlabel('Subject'); ylabel('Throw time (seconds)'); ylim([0.8, 3.5]);
set(findobj(gca,'type','line'),'linew',.75)

tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% plot single-trial distance X condition
clc; close all;
r = struct2table(load_r.r) ;
fig = figure('Position', [0 0 1366 786], 'Visible', 'off'); 
h1 = histogram(r.distance(~~r.pres,:), 'BinEdges',1:1:49, 'FaceAlpha',0.15); hold on;
h2 = histogram(r.distance(~r.pres,:), 'BinEdges',1:1:49, 'FaceAlpha',0.15);
h1v = h1.Values;
h2v = h2.Values; clf;
bar(2:1:49,[h1v;h2v]', 1,'stacked', 'FaceAlpha',0.15);
ylabel('Trials');
xlabel('Distance (sqrt(cm))');
xline(mean(r.distance(~~r.pres,:)), 'LineWidth',2, 'Color','b');
xline(mean(r.distance(~r.pres,:)), 'LineWidth',2, 'Color','r');
set(gca,'FontSize',16);
t = title('Dart''s Distance from Bull''s-Eye by Condition');
legend({'Target Present','Target Absent', 'Target Present Mean', 'Target Absent Mean'});
tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% plot single-trial throwtime X condition with subject means
clc; close all;
r = struct2table(load_r.r) ;
figure('Position', [0 0 1024 786]);
% figure('Position', [0 0 1024 786], 'Visible','off');
h1 = histogram(r.throwtime(~~r.pres,:), 'BinEdges',0.9:.05:1.8, 'FaceAlpha',0.15); hold on;
h2 = histogram(r.throwtime(~r.pres,:), 'BinEdges',0.9:.05:1.8, 'FaceAlpha',0.15);
h1v = h1.Values;
h2v = h2.Values; clf;
bar(.95:.05:1.8,[h1v;h2v]', 1,'stacked', 'FaceAlpha',0.15);
ylabel('Trials');
xlabel('Throw time (sqrt(seconds))');
xline(mean(r.throwtime(~~r.pres,:)), 'LineWidth',2, 'Color','b');
xline(mean(r.throwtime(~r.pres,:)), 'LineWidth',2, 'Color','r');
set(gca,'FontSize',16);
t = title('Time Taken to Throw by Condition');
legend({'Target Present','Target Absent', 'Target Present Mean', 'Target Absent Mean'});
tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% distance vs condition statistics
% single-trial
r = struct2table(load_r.r) ;
ta = r.distance(r.pres == 0);
tp = r.distance(r.pres == 1);
[h,p] = perm_ttest2(ta,tp,'nperm', 10000)

% summary statistic
abs_subj_means = []; pres_subj_means = [];
for subj_i = 1:length(subjs_to_include)
	abs_subj_means(end+1,:) = mean(r.distance(r.subject == subj_i & r.pres == 1));
	pres_subj_means(end+1,:) = mean(r.distance(r.subject == subj_i & r.pres == 0));
end
[h,p] = perm_ttest2(abs_subj_means,pres_subj_means,'nperm', 10000)

%% throwtime vs condition statistics
% single-trial
r = struct2table(load_r.r) ;
ta = r.throwtime(r.pres == 0);
tp = r.throwtime(r.pres == 1);
[h,p] = perm_ttest2(ta,tp,'nperm', 10000)

% summary statistic
abs_subj_means = []; pres_subj_means = [];
for subj_i = 1:length(subjs_to_include)
	abs_subj_means(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.pres == 1));
	pres_subj_means(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.pres == 0));
end
[h,p] = perm_ttest2(abs_subj_means,pres_subj_means,'nperm', 10000)

%% plot single-trial distance vs delay by condition
clc; fig = figure(3); clf;
set(fig, 'Position', [0 0 1024 786]*2);

r = struct2table(load_r.r) ;
abs_delay = r.delay(r.pres==0);
pres_delay = r.delay(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);

% single-trial statistics
disp('Single Trial')
disp('Target Absent')
[rho, p] = mperm_corr(abs_delay,abs_distance)
disp('Target Present')
[rho, p] = mperm_corr(pres_delay,pres_distance)

% single trial plot
subplot(2,1,1)
plot(abs_delay, abs_distance,'.'); hold on;
title('Target Absent');
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay-time (seconds)')
ylabel('Distance (sqrt(cm))')
lsline;
set(gca,'FontSize',11)

subplot(2,1,2)
plot(pres_delay, pres_distance,'.'); hold on;
title('Target Present');
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay-time (seconds)')
ylabel('Distance (sqrt(cm))')
lsline;
set(gca,'FontSize',11)

%% summary statistic correlation using subject means
abs_subj_means = []; abs_subj_delays = [];
pres_subj_means = []; pres_subj_delays = [];
for subj_i = 1:length(subjs_to_include)
	for delay = 3:9
		abs_subj_delays(end+1,:) = delay;
		abs_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		rm_abs_subj_means(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		pres_subj_delays(end+1,:) = delay;
		pres_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==1));
		rm_abs_subj_means(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
	end
end
disp('Summary Statistics')
disp('Target Absent')
[rho,p] = mperm_corr(pres_subj_delays,pres_subj_means)
disp('Target Present')
[rho,p] = mperm_corr(abs_subj_delays,abs_subj_means)

% summary statistic plot using subject means
subplot(2,1,1)
plot(abs_subj_delays,abs_subj_means,'r*')
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay-time (seconds)')
ylabel('Distance (sqrt(cm))')
lsline;
set(gca,'FontSize',11)

subplot(2,1,2)
plot(pres_subj_delays,pres_subj_means,'r*')
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay-time (seconds)')
ylabel('Distance (sqrt(cm))')
lsline;
set(gca,'FontSize',11)

% subplot(3,1,3)
% plot(pres_subj_delays,abs_subj_means-pres_subj_means,'r*')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
% xlabel('Delay-time (seconds)')
% ylabel('Distance (sqrt(cm))')
% lsline;
% set(gca,'FontSize',11)

% summary statistic correlation using subject x cond means
abs_cond_subj_means = [];
pres_cond_subj_means = [];
for delay = 3:9
	abs_cond_subj_means(end+1,:) = mean(abs_subj_means(abs_subj_delays==delay));
	pres_cond_subj_means(end+1,:) = mean(pres_subj_means(pres_subj_delays==delay));
end
disp('Summary Statistics')
disp('Target Absent')
[rho,p] = mperm_corr((3:9)',abs_cond_subj_means)
disp('Target Present')
[rho,p] = mperm_corr((3:9)',pres_cond_subj_means)
disp('Target Absent - Target Present')
[rho,p] = mperm_corr((3:9)',abs_cond_subj_means-pres_cond_subj_means)

% summary statistic plot using subject x cond means
subplot(2,1,1)
plot(3:9,abs_cond_subj_means,'ko')
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay-time (seconds)')
ylabel('Distance (sqrt(cm))')
lsline;
set(gca,'FontSize',11)

subplot(2,1,2)
plot(3:9,pres_cond_subj_means,'ko')
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Delay-time (seconds)')
ylabel('Distance (sqrt(cm))')
lsline;
set(gca,'FontSize',11)

%%
% fitlme(r,'distance ~ throwtime*delay*pres')
r.pres = categorical(r.pres);
r.delay = r.delay-mean(unique(r.delay));
lme = fitlme(r,'distance ~ throwtime*delay*pres + (1+delay|subject)+ (1+throwtime|subject)', 'FitMethod', 'REML', 'StartMethod','random')

% fitlme(r,'distance ~ (throwtime*pres|subject) + (throwtime*pres*delay|subject)')

%% distance vs throwtime
clc; fig = figure(4); clf;
set(fig, 'Position', [0 0 1024 786]*2);

for mem = 0:1
	r = struct2table(load_r.r) ;
	r(r.pres==mem,:)=[];
	
	% single trial plot and correlation
	subplot(2,1,mem+1);
	if mem == 1; disp('Target Absent');  title('Target Absent'); end; hold on;
	if mem == 0; disp('Target Present'); title('Target Present'); end; hold on;
	plot(r.distance, r.throwtime,'.'); hold on;
	xlim([0, 8]); 	ylim([0.75, 2])
	ylabel('Distance (sqrt(cm))')
	xlabel('Throw-time (sqrt(seconds))')
	disp('Single Trial')
	[rho, p] = mperm_corr(r.throwtime,r.distance);
	[rho, p] 

	% summary statistic plot and correlation
	subplot(2,1,mem+1);
	for subj_i = 1:length(subjs_to_include)
		subj_is = r.subject==subj_i;
		mean_throwtime = mean(r.throwtime(subj_is));
		mean_distance = mean(r.distance(subj_is));
		means(subj_i,:) = [mean_throwtime, mean_distance];
	end
	plot(means(:,2),means(:,1),'r*'); hold on;
	% set axes and least squares line
	lsline;
	xlim([0, 8]); 	ylim([0.75, 2])
	set(gca,'FontSize',11)
	disp('Summary Statistics')
	[rho, p] =mperm_corr(means(:,1),means(:,2));
	[rho, p]
	
end

%% make regression vars
clc;
r = struct2table(load_r.r) ;
r(any(isnan(r.Variables),2),:) = [] ; % remove all nans

% % within-subject normalize distance
% distance = r.distance.^(1/2) ;% sqrt(paper units)
% % distance = r.distance ;% sqrt(paper units)
% for subj_i = 1:length(subjs_to_include)
% 	subj_id = subjs_to_include(subj_i) ;
% 	subj_inds = r.subject==str2double(subj_id{:}) ;
% 	r.distance(subj_inds) =  normalize(r.distance(subj_inds)) ;
% end

% normalize the rest and clean
pres = r.pres ;
subject = r.subject  ;
r.pres = [] ;
r.subject = [];
r.Variables = normalize(r.Variables) ;
rej_i = any(abs(r.Variables) > 6,2) ; % set rejection criterion here
% rej_i = false(height(r)) ; % set rejection criterion here
r(rej_i,:) = [] ; sum(rej_i);
r.trialnum = [] ; r.delay = [] ; r.delaystilltime = [] ;
r.octant = [] ;  r.position = [] ;
r.delay_Front_gamma = [] ; r.delay_Back_gamma = [] ;
r.pre_throw_Front_gamma = [] ; r.pre_throw_Back_gamma = [] ;
r.pres = categorical(pres(~rej_i)) ;
r.subject = categorical(subject(~rej_i));

%% create variables lists
chan_labs = {'Front','Back'} ;
freq_labs = {'theta', 'alpha'} ;
dvarnames = {} ; pvarnames = {} ;
for chan_lab = chan_labs
	for freq_lab = freq_labs
		pvarnames{end+1} = ['pre_throw_',chan_lab{:},'_',freq_lab{:}] ;
		dvarnames{end+1} = ['delay_',chan_lab{:},'_',freq_lab{:}] ;
	end
end
bvarnames = {'pres'} ;

%% create time-averaged eeg columns
for field = pvarnames
	r.(field{:}) = normalize(r.(field{:})./r.throwtime) ;
end
for field = dvarnames
	r.(field{:}) = normalize(r.(field{:})./r.delaytime) ;
end
r.throwtime = [] ;
r.delaytime = [] ;
lin_r = r ;
int_r = table2struct(lin_r) ;
int2_r = table2struct(lin_r) ;

%% create coupling terms
vars = {pvarnames{:}}' ;
vcs = nchoosek(vars,2) ;
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1} ;
	v2n = vcs{vc_i,2} ;
	v1 = num2cell(([r.pres].*[r.(v1n)])) ;
	v2 = num2cell(([r.pres].*[r.(v2n)])) ;
	v12 = num2cell(([r.pres].*[r.(v1n)].*[r.(v2n)])) ;
	[int_r.(['presx',v1n])] = v1{:} ;
	[int_r.(['presx',v2n])]= v2{:} ;
	[int2_r.(['presx',v1n])] = v1{:} ;
	[int2_r.(['presx',v2n])]= v2{:} ;
	[int2_r.(['presx',v1n,'x',v2n])] = v12{:} ;
end
vars = {dvarnames{:}}' ;
vcs = nchoosek(vars,2) ;
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1} ;
	v2n = vcs{vc_i,2} ;
	v1 = num2cell([r.pres].*[r.(v1n)]) ;
	v2 = num2cell([r.pres].*[r.(v2n)]) ;
	v12 = num2cell([r.pres].*[r.(v1n)].*[r.(v2n)]) ;
	[int_r.(['presx',v1n])] = v1{:} ;
	[int_r.(['presx',v2n])]= v2{:} ;
	[int2_r.(['presx',v1n])] = v1{:} ;
	[int2_r.(['presx',v2n])]= v2{:} ;
	[int2_r.(['presx',v1n,'x',v2n])] = v12{:} ;
end
int_r = struct2table(int_r) ;
int2_r = struct2table(int2_r) ;

%% get reference diagnostic measures for regression
y = [r.distance]' ;
ref = [sum((y-mean(y)).^2),...
	sum(abs(y-mean(y))),...
	]' ;
ref_mean_press = 0 ;
for i = 1:length(y)
	pop_y = y ;
	pop_y(i) = [] ;
	ref_mean_press = ref_mean_press + (y(i)-mean(pop_y)).^2 ;
end
ref_median_press = 0 ;
for i = 1:length(y)
	pop_y = y ;
	pop_y(i) = [] ;
	ref_median_press = ref_median_press + (y(i)-median(pop_y)).^2 ;
end

%% regression
% close all ;
% figure(1) ; figure(2) ;
% pause(1)
% clc
int = false ; % should be false with pres as a dummy variable
% r = table2struct(int_r) ; % linear+int model
% r = table2struct(int2_r) ; % linear+second order int model


%% linear mixed effects
tab_r = lin_r;
tab_r.pres = logical(double(tab_r.pres))

% tab_r.abs = ~logical(double(tab_r.pres))
% mdl = fitlmematrix(X,y,'VarNames',fields, 'ResponseVar','distance', 'intercept', int) ;
% form = {
% 	'distance ~ 1 +'
% 	'pres*delay_Front_theta +'
% 	'pres*pre_throw_Front_theta +'
% 	'pres*delay_Front_alpha +'
% 	'pres*pre_throw_Front_alpha +'
% 	'pres*delay_Back_theta +'
% 	'pres*pre_throw_Back_theta +'
% 	'pres*delay_Back_alpha +'
% 	'pres*pre_throw_Back_alpha +'
% 	'(pres*delay_Front_theta|subject) +'
% 	'(pres*pre_throw_Front_theta|subject)  +'
% 	'(pres*delay_Front_alpha|subject)  +'
% 	'(pres*pre_throw_Front_alpha|subject)  +'
% 	'(pres*delay_Back_theta|subject)  +'
% 	'(pres*pre_throw_Back_theta|subject)  +'
% 	'(pres*delay_Back_alpha|subject)  +'
% 	'(pres*pre_throw_Back_alpha|subject)'
% 	};
% form = {
% 	'distance ~ '
% 	'pres*delay_Front_theta +'
% 	'pres*pre_throw_Front_theta +'
% 	'pres*delay_Front_alpha +'
% 	'pres*pre_throw_Front_alpha +'
% 	'pres*delay_Back_theta +'
% 	'pres*pre_throw_Back_theta +'
% 	'pres*delay_Back_alpha +'
% 	'pres*pre_throw_Back_alpha +'
% 	'(1+pres|subject)'
% 	};
% form = {
% 	'distance ~ '
% 	'pres*delay_Front_theta +'
% 	'pres*pre_throw_Front_theta +'
% 	'pres*delay_Back_alpha +'
% 	'pres*pre_throw_Back_alpha +'
% 	'(1+pres|subject)'
% 	};
form = {
	'delay_Front_theta ~'
	'pres'
	};
mdl = fitlme(tab_r, ...
	[form{:}], ...
	'FitMethod', 'REML',...
	'DummyVarCoding','full',...
	'StartMethod','random');
[[sum((mdl.Residuals.Raw).^2),...
	sum(abs(mdl.Residuals.Raw))]',ref]
mdl

%%
clear
ds = dataset();
% ds.y = rand(100,1);
res = randn(100,1);
ds.X = ones(100,1)+res;
ds.A = randi([0 1],100,1);
eps = randn(100,1);
b = .2*randn(100,1);
ds.y = 2.*ds.X + .3.*eps.*ds.A+ .7.*eps.*~ds.A;
fitlme(ds,'y~X+(X|A)')
%% some random data.
ds.y = rand(8,1);
ds.A = [true;true;true;true;false;false;false;false];
ds.B = rand(8,1);
ds.C = categorical([true;true;true;true;false;false;false;false]);
lme = fitlme(ds,'y ~ 1 + B + C')

%%
tab_r = lin_r;
y = tab_r.delay_Front_theta;
X = [double(tab_r.pres)];
Z = [ones(height(tab_r),1)];
G = tab_r.subject;
mdl = fitlmematrix(X,y,Z, G)
%% non-robust
%%
tab_r = lin_r;
tab_r.subject = [];
tab_r.pres = double(tab_r.pres);
r = table2struct(tab_r) ; % linear model
y = [r.distance]' ;
X = r ;
fields = fieldnames(X) ;
X = rmfield(X,'distance') ;
X = cell2mat(struct2cell(X))' ;
mdl = fitlm(X,y,'VarNames',fields, 'ResponseVar','distance', 'intercept', int) ;
[[sum((mdl.Residuals.Raw).^2),...
	sum(abs(mdl.Residuals.Raw))]',ref]
ref_mean_press
press_var = press([X,y])
mdl
%% non-robust bootstrap
mdl = fitlm(X,y,'VarNames',fields, 'ResponseVar','distance', 'intercept', int) ;
bhs = runboot(X,y,mdl.Residuals.Raw,100000) ;
quantile(bhs,[0.05, .95]) ;
figure(1) ;
ncoef = size(bhs,2) ;
for i = 1:ncoef
	subplot( round(sqrt(ncoef)) , ceil(sqrt(ncoef)) , i) ; hist(bhs(:,i))
end

% robust
mdl = fitlm(X,y,'VarNames',fields, 'RobustOpts','on', 'ResponseVar','distance', 'intercept', int) ;
[sum((mdl.Residuals.Raw).^2)<ref,...
	sum(abs(mdl.Residuals.Raw))<ref]'
press_var = 0 ;
for i = 1:length(y)
	if mod(i,100) == 0
		disp([num2str(100*(i/length(y))),'% done']) ;
	end
	ind = true(size(y)) ;
	ind(i) = false ;
	mdl = fitlm(X(ind,:),y(ind),'VarNames',fields, 'RobustOpts','on', 'ResponseVar','distance', 'intercept',int) ;
	rb =  mdl.Coefficients.Estimate ;
	press_var = press_var + ((y(i) - X(i,:)*rb).^2) ;
end
ref_mean_press
press_var

% robust bootstrap
mdl = fitlm(X,y,'VarNames',fields, 'RobustOpts','on', 'ResponseVar','distance', 'intercept', int) ;
bhs = runboot(X,y,mdl.Residuals.Raw,100000) ;
quantile(bhs,[0.05, .95]) ;
figure(2) ;
ncoef = size(bhs,2) ;
for i = 1:ncoef
	subplot( round(sqrt(ncoef)) , ceil(sqrt(ncoef)) , i) ; hist(bhs(:,i))
end