%% init
clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
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
	};
srate = 512;

% make regression vars
load('ic_r.mat');
r = struct2table(r);
r(any(isnan(r.Variables),2),:) = [];

% within-subject normalize distance
distance = r.distance.^(1/2);% sqrt(paper units)
% distance = r.distance;% sqrt(paper units)
for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include(subj_i);
	subj_inds = r.subject==str2double(subj_id{:});
	r.distance(subj_inds) =  normalize(r.distance(subj_inds));
end

% normalize the rest and clean
pres = r.pres;
r.pres = [];
r.Variables = normalize(r.Variables);
% rej_i = any(abs(r.Variables)>10,2);
% sum(rej_i)
rej_i = false(size(pres));
r(rej_i,:) = [];
r.trialnum = []; r.delay = []; r.delaystilltime = [];
r.octant = []; r.subject = []; r.position = [];
r.delay_Front_gamma = []; r.delay_Back_gamma = [];
r.pre_throw_Front_gamma = []; r.pre_throw_Back_gamma = [];

r.pres = pres(~rej_i);

% create variables lists
chan_labs = {'Front','Back'};
freq_labs = {'theta', 'alpha'};
dvarnames = {}; pvarnames = {};
for chan_lab = chan_labs
	for freq_lab = freq_labs
		pvarnames{end+1} = ['pre_throw_',chan_lab{:},'_',freq_lab{:}];
		dvarnames{end+1} = ['delay_',chan_lab{:},'_',freq_lab{:}];
	end
end
bvarnames = {'pres'};

% create time-averaged eeg columns
for field = pvarnames
	r.(field{:}) = r.(field{:})./r.throwtime;
end
for field = dvarnames
	r.(field{:}) = r.(field{:})./r.delaytime;
end
r.throwtime = [];
r.delaytime = [];
lin_r = r;
r = table2struct(r);
% create coupling terms
vars = {pvarnames{:}}';
vcs = nchoosek(vars,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	v1 = num2cell([r.pres].*[r.(v1n)]);
	v2 = num2cell([r.pres].*[r.(v2n)]);
	v12 = num2cell([r.pres].*[r.(v1n)].*[r.(v2n)]);
	[r.(['presx',v1n])] = v1{:};
	[r.(['presx',v2n])]= v2{:};
	[r.(['presx',v1n,'x',v2n])] = v12{:};
end
vars = {dvarnames{:}}';
vcs = nchoosek(vars,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	v1 = num2cell([r.pres].*[r.(v1n)]);
	v2 = num2cell([r.pres].*[r.(v2n)]);
	v12 = num2cell([r.pres].*[r.(v1n)].*[r.(v2n)]);
	[r.(['presx',v1n])] = v1{:};
	[r.(['presx',v2n])]= v2{:};
	[r.(['presx',v1n,'x',v2n])] = v12{:};
end
r = struct2table(r);
int_r = r;

clc
y = [r.distance]';
ref = [sum((y-mean(y)).^2),...
	sum((y-median(y)).^2),...
	sum(abs(y-mean(y))),...
	sum(abs(y-median(y)))...
	]';
ref_press = 0;
for i = 1:length(y)
	pop_y = y;
	pop_y(i) = [];
	ref_press = ref_press + (y(i)-mean(pop_y)).^2;
end
% mdl = fitlm(r,'distance~1');
% r_median_absolute_deviation = nanmedian(abs(mdl.Residuals.Raw));
% r_mean_absolute_deviation = nanmean(abs(mdl.Residuals.Raw));
% r_root_mean_squared_error = mdl.RMSE;
% r_root_median_squared_error = sqrt(nanmedian(mdl.Residuals.Raw.^2));
% ref = [r_median_absolute_deviation;
% 	r_mean_absolute_deviation;
% 	r_root_mean_squared_error;
% 	r_root_median_squared_error]
% y = [r.distance]';
%
% mdl = fitlm(ones(size(y))*median(y),y,'Intercept',false,'RobustOpts','on');
% r_median_absolute_deviation = nanmedian(abs(mdl.Residuals.Raw));
% r_mean_absolute_deviation = nanmean(abs(mdl.Residuals.Raw));
% r_root_mean_squared_error = mdl.RMSE;
% r_root_median_squared_error = sqrt(nanmedian(mdl.Residuals.Raw.^2));
% rb_ref = [r_median_absolute_deviation;
% 	r_mean_absolute_deviation;
% 	r_root_mean_squared_error;
% 	r_root_median_squared_error]
%% linear model
clc
r = table2struct(lin_r);
y = [r.distance]';
X = r;
fields = fieldnames(X);
X = rmfield(X,'distance');
X = cell2mat(struct2cell(X))';

% non-robust
mdl = fitlm(X,y,'VarNames',fields, 'ResponseVar','distance');
ref
[sum((mdl.Residuals.Raw).^2)<ref,...
	sum(abs(mdl.Residuals.Raw))<ref]'
press([X,y])<ref_press

%% robust
mdl = fitlm(X(ind,:),y(ind),'VarNames',fields, 'RobustOpts','on', 'ResponseVar','distance', 'intercept',false);

press_var = 0;
for i = 1:length(y)
	if mod(i,100) == 0; disp([num2str(100*(i/length(y))),'% done']);	end
	ind =true(size(y));
	ind(i) = false;
	mdl = fitlm(X(ind,:),y(ind),'VarNames',fields, 'RobustOpts','on', 'ResponseVar','distance', 'intercept',false);
	rb =  mdl.Coefficients.Estimate;
	press_var = press_var + (y(i) - X(i,:)*rb);
end
press_var<ref_press
[sum((mdl.Residuals.Raw).^2)<ref,...
	sum(abs(mdl.Residuals.Raw))<ref]'
%% linear+interaction model
clc
r = table2struct(int_r);
y = [r.distance]';
X = r;
fields = fieldnames(X);
X = rmfield(X,'distance');
X = cell2mat(struct2cell(X))';
b = X\y;
res = X*b;
mdl = fitlm(X,y,'VarNames',fields, 'ResponseVar','distance');
[r_mad, r_rmse]
[nanmedian(abs(mdl.Residuals.Raw)),...
	mdl.RMSE,...
	mdl.Rsquared.Adjusted,...
	coefTest(mdl)...
	]
mdl = fitlm(X,y,'VarNames',fields, 'RobustOpts','on', 'ResponseVar','distance');
[rb_mad, rb_rmse]
[nanmedian(abs(mdl.Residuals.Raw)),...
	mdl.RMSE,...
	mdl.Rsquared.Adjusted,...
	coefTest(mdl)...
	]

%% linear+interaction model
clc
r = table2struct(lin_r);
y = [r.distance]';
X = r;
fields = fieldnames(X);
X = rmfield(X,'distance');
X = cell2mat(struct2cell(X))';
b = X\y;
res = X*b;
mdl = fitlm(X,y,'VarNames',fields);
[r_mad, r_rmse]
[nanmedian(abs(mdl.Residuals.Raw)),...
	mdl.RMSE,...
	mdl.Rsquared.Adjusted,...
	coefTest(mdl)...
	]
mdl = fitlm(X,y,'VarNames',fields, 'RobustOpts','on');
[rb_mad, rb_rmse]
[nanmedian(abs(mdl.Residuals.Raw)),...
	mdl.RMSE,...
	mdl.Rsquared.Adjusted,...
	coefTest(mdl)...
	]

%% linear+interaction model
clc
r = lin_r;
vcn = {bvarnames{:}, pvarnames{:}, dvarnames{:}};
vars = {pvarnames{:}}';
vcs = nchoosek(vars,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	vcn{end+1} = ['pres:',v1n];
	vcn{end+1} = ['pres:',v2n];
	vcn{end+1} = ['pres:',v1n,':',v2n];
end

vars = {dvarnames{:}}';
vcs = nchoosek(vars,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	vcn{end+1} = ['pres:',v1n];
	vcn{end+1} = ['pres:',v2n];
	vcn{end+1} = ['pres:',v1n,':',v2n];
end

eqn = 'distance ~ ';
eqn = [eqn, strjoin(vcn,' + ')];

mdl = fitlm(r,eqn);
[r_mad, r_rmse]
[nanmedian(abs(mdl.Residuals.Raw)),...
	mdl.RMSE,...
	mdl.Rsquared.Adjusted,...
	coefTest(mdl)...
	]
mdl = fitlm(r,eqn, 'RobustOpts','on');
x = (mdl.Formula.Terms*table2array(mdl.Variables)')';
[~,c,~] = collintest(x,'Display','off');
[rb_mad, rb_rmse]
[nanmedian(abs(mdl.Residuals.Raw)),...
	mdl.RMSE,...
	mdl.Rsquared.Adjusted,...
	coefTest(mdl)...
	]

%% linear+interaction bootstrap
clc
x = (mdl.Formula.Terms*table2array(mdl.Variables)')';
bhs = runboot(x,r.distance,mdl.Residuals.Raw,10000);
