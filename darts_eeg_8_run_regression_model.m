%% init
clear ; close all ; clc ;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/' ;
cd(script_dir) ;
addpath([script_dir,'deps/']) ;
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

%% make regression vars
clc;
r = struct2table(load_r.r) ;
r(any(isnan(r.Variables),2),:) = [] ;

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
	sum((y-median(y)).^2),...
	sum(abs(y-mean(y))),...
	sum(abs(y-median(y)))...
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
tab_r = lin_r;
% r = table2struct(lin_r) ; % linear model
% y = [r.distance]' ;
% X = r ;
% fields = fieldnames(X) ;
% X = rmfield(X,'distance') ;
% X = cell2mat(struct2cell(X))' ;

%% mixed effects
% mdl = fitlmematrix(X,y,'VarNames',fields, 'ResponseVar','distance', 'intercept', int) ;
mdl = fitlme(tab_r, 'distance ~ 1+pres + delay_Front_theta+(delay_Front_theta|subject)', 'FitMethod', 'REML')
mdl.Rsquared.Adjusted
%% non-robust
mdl = fitlm(X,y,'VarNames',fields, 'ResponseVar','distance', 'intercept', int) ;
[sum((mdl.Residuals.Raw).^2)<ref,...
	sum(abs(mdl.Residuals.Raw))<ref]'
ref_mean_press
press_var = press([X,y])

% non-robust bootstrap
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