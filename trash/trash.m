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

%% pca on linear eeg terms
r = table2struct(lin_r);
fields = {dvarnames{:},pvarnames{:}};
[t.(fields{:})] = [r.(fields{:})];

%%
% eval( ['leeg_table = table(', strjoin([dvarnames, pvarnames],','),');'] );
[lcoeff,lscore,~,~,lexplained,~] = pca(table2array(leeg_table));
eval( ['r = table(', strjoin([bvarnames],','),');'] );

r.lpc1= lscore(:,1);
r.lpc2 = lscore(:,2);
r.lpc3 = lscore(:,3);

% robust, no interactions
mdl_robust = fitlm(r,'ResponseVar','distance', 'RobustOpts','on')
mad_robust = nanmedian(abs(mdl_robust.Residuals.Raw))

% non-robust, no interactions
mdl = fitlm(r,'ResponseVar','distance')
mad = nanmedian(abs(mdl.Residuals.Raw))

%% pca on just coupling terms
eval( ['deeg_table = table(', strjoin(dvarnames,','),');'] );
eval( ['peeg_table = table(', strjoin(pvarnames,','),');'] );
t = deeg_table;
dpca = table();
vcs = nchoosek(t.Properties.VariableNames,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	vcn = [v1n,'x',v2n];
	evalc(['dpca.',vcn,' = t.',v1n,'.*t.',v2n,';']);
end

t = peeg_table;
ppca = table();
vcs = nchoosek(t.Properties.VariableNames,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	vcn = [v1n,'x',v2n];
	evalc(['ppca.',vcn,' = t.',v1n,'.*t.',v2n,';']);
end

[dcoeff,dscore,~,~,dexplained,~] = pca(table2array(dpca));
[pcoeff,pscore,~,~,pexplained,~] = pca(table2array(ppca));
eval( ['r = table(', strjoin([bvarnames,dvarnames, pvarnames],','),');'] );
r.dpc1= dscore(:,1);
r.dpc2 = dscore(:,2);
r.ppc1= pscore(:,1);
r.ppc2 =pscore(:,2);

eqn = [...
	'distance ~',...
	'pres +',...
	'throwtime +',...
	'delaytime +',...
	'pres:throwtime +',...
	'pres:delaytime +',...
	'throwtime:delaytime +',...
	'pres:throwtime:delaytime +',...
	'dpc1 +',...
	'dpc2 +',...
	'dpc1:dpc2 +',...
	'ppc1 +',...
	'ppc2 +',...
	'ppc1:ppc2 +',...
	'pres:dpc1 +',...
	'pres:dpc2 +',...
	'pres:dpc1:dpc2 +',...
	'pres:ppc1 +',...
	'pres:ppc2 +',...
	'pres:ppc1:ppc2 +',...
	'delaytime:dpc1 +',...
	'delaytime:dpc2 +',...
	'delaytime:dpc1:dpc2 +',...
	'throwtime:ppc1 +',...
	'throwtime:ppc2 +',...
	'throwtime:ppc1:ppc2 +',...
	'pres:delaytime:dpc1 +',...
	'pres:delaytime:dpc2 +',...
	'pres:delaytime:dpc1:dpc2 +',...
	'pres:throwtime:ppc1 +',...
	'pres:throwtime:ppc2 +',...
	'pres:throwtime:ppc1:ppc2'...
	];

% robust eqn
mdl_robust = fitlm(r,eqn, 'RobustOpts','on')
mad_robust = nanmedian(abs(mdl_robust.Residuals.Raw))

% non-robust eqn
% mdl = fitlm(r,eqn)
% mad = nanmedian(abs(mdl.Residuals.Raw))

%% pca on linear and coupling terms
eval( ['deeg_table = table(', strjoin(dvarnames,','),');'] );
eval( ['peeg_table = table(', strjoin(pvarnames,','),');'] );
eval( ['aeeg_table = table(', strjoin([dvarnames,pvarnames],','),');'] );

t = deeg_table;
vcs = nchoosek(t.Properties.VariableNames,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	vcn = [v1n,'x',v2n];
	evalc(['deeg_table.',vcn,' = t.',v1n,'.*t.',v2n,';']);
end

t = peeg_table;
vcs = nchoosek(t.Properties.VariableNames,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	vcn = [v1n,'x',v2n];
	evalc(['peeg_table.',vcn,' = t.',v1n,'.*t.',v2n,';']);
end

t = aeeg_table;
vcs = nchoosek(t.Properties.VariableNames,2);
for vc_i = 1:length(vcs)
	v1n = vcs{vc_i,1};
	v2n = vcs{vc_i,2};
	vcn = [v1n,'x',v2n];
	evalc(['aeeg_table.',vcn,' = t.',v1n,'.*t.',v2n,';']);
end

[dcoeff,dscore,~,~,dexplained,~] = pca(table2array(deeg_table));
[pcoeff,pscore,~,~,pexplained,~] = pca(table2array(peeg_table));
[acoeff,ascore,~,~,aexplained,~] = pca(table2array(aeeg_table));

dpc1= dscore(:,1);
dpc2 = dscore(:,2);
ppc1= pscore(:,1);
ppc2 =pscore(:,2);
apc1= ascore(:,1);
apc2 = ascore(:,2);
apc3 = ascore(:,3);

% without cross-temporal interactions
varnames = {'distance'; 'throwtime'; 'delaytime';'pres'};
varnames{end+1} = 'dpc1';
varnames{end+1} = 'dpc2';
varnames{end+1} = 'ppc1';
varnames{end+1} = 'ppc2';
eval( ['r = table(', strjoin(varnames,','),');'] );

% with cross-temporal interactions
varnames = {'distance'; 'throwtime'; 'delaytime';'pres'};
varnames{end+1} = 'apc1';
varnames{end+1} = 'apc2';
varnames{end+1} = 'apc3';
eval( ['r = table(', strjoin(varnames,','),');'] );

% run as non-eqn regression

%% reduced pca eqn
eqn = [...
	'distance ~',...
	'pres +',...
	'throwtime +',...
	'delaytime +',...
	'pres:throwtime +',...
	'pres:delaytime +',...
	'throwtime:delaytime +',...
	'pres:throwtime:delaytime +',...
	'pres:dpc1 +',...
	'pres:dpc2 +',...
	'pres:ppc1 +',...
	'pres:ppc2 +',...
	'delaytime:dpc1 +',...
	'delaytime:dpc2 +',...
	'throwtime:ppc1 +',...
	'throwtime:ppc2 +',...
	'pres:delaytime:dpc1 +',...
	'pres:delaytime:dpc2 +',...
	'pres:throwtime:ppc1 +',...
	'pres:throwtime:ppc2'];

%% reduced pca eqn, cross-temporal
eqn = [...
	'distance ~',...
	'pres +',...
	'throwtime +',...
	'delaytime +',...
	'pres:throwtime +',...
	'pres:delaytime +',...
	'throwtime:delaytime +',...
	'pres:throwtime:delaytime +',...
	'pres:dpc1 +',...
	'pres:dpc2 +',...
	'pres:ppc1 +',...
	'pres:ppc2 +',...
	'ppc1:dpc1 +',...
	'ppc1:dpc2 +',...
	'ppc2:dpc1 +',...
	'ppc2:dpc2 +',...
	'delaytime:dpc1 +',...
	'delaytime:dpc2 +',...
	'throwtime:ppc1 +',...
	'throwtime:ppc2 +',...
	'pres:delaytime:dpc1 +',...
	'pres:delaytime:dpc2 +',...
	'pres:throwtime:ppc1 +',...
	'pres:throwtime:ppc2'];

%% behavioral eqn
eqn = [...
	'distance ~',...
	'pres +',...
	'throwtime +',...
	'delaytime +',...
	'pres:throwtime +',...
	'pres:delaytime +',...
	'throwtime:delaytime +',...
	'pres:throwtime:delaytime'];

%% robust, no interactions
mdl_robust = fitlm(r,'ResponseVar','distance', 'RobustOpts','on')
mad_robust = nanmedian(abs(mdl_robust.Residuals.Raw))

%% non-robust, no interactions
mdl = fitlm(r,'ResponseVar','distance')
mad = nanmedian(abs(mdl.Residuals.Raw))

%% robust eqn
mdl_robust = fitlm(r,eqn, 'RobustOpts','on')
mad_robust = nanmedian(abs(mdl_robust.Residuals.Raw))

%% non-robust eqn
mdl = fitlm(r,eqn)
mad = nanmedian(abs(mdl.Residuals.Raw))

%% robust interactions
mdl_robust = fitlm(r,'interactions','ResponseVar', 'distance', 'RobustOpts','on')
mad_robust = nanmedian(abs(mdl_robust.Residuals.Raw))

%% non-robust interactions
mdl = fitlm(r,'interactions','ResponseVar', 'distance')
mad = nanmedian(abs(mdl.Residuals.Raw))

%% predict condition
bvarnames = {'pres'};
eval( ['r = table(', strjoin([bvarnames,dvarnames],','),');'] );
mdl = fitglm(r,'interactions','ResponseVar', 'pres')
% mdl = stepwiseglm(r,'interactions','ResponseVar', 'pres')
%
% mdl.step()
% mdl = fitglm(r,'pres~delay_Fz_theta+delay_Cz_theta')

%% can eeg explain distance
bvarnames = {'distance'};
eval( ['r = table(', strjoin([bvarnames,dvarnames,pvarnames],','),');'] );
% mdl = fitglm(r,'distance~delay_Fz_theta*delay_Cz_theta')
mdl = stepwiseglm(r,'interactions','ResponseVar', 'distance', 'Criterion','AdjRsquared')
mad = nanmedian(abs(mdl.Residuals.Raw))

%% univariate
% mdl = fitglm(r,'distance~pres')
distcoup = fitglm(r,'distance~delay_Fz_gamma:delay_Cz_alpha')

% mdl = fitglm(r,'pres~delay_Fz_theta+delay_Cz_theta+delay_Pz_theta+delay_Oz_alpha')

% mdl = fitlm(r,'distance~pres+pres:delaytime:delay_Fz_theta+pres:throwtime:pre_throw_Oz_alpha+pres:throwtime:pre_throw_Cz_gamma')
% mdl = fitlm(r,'distance~pres:delaytime:delay_Cz_alpha')
% mdl = fitlm(r,'distance~delaytime:delay_Fz_theta+delaytime:delay_Cz_theta+delaytime:delay_Pz_theta+delaytime:delay_Oz_theta')

% mdl = fitlm(r,'distance~pres:throwtime:pre_throw_Oz_alpha')

%% clean results
summary_results = sortrows(anova(mdl_robust),5,'ascend');
coeff_results = sortrows(mdl_robust.Coefficients,4,'ascend');
save('results_unz.mat');

% save results to xls file for table and fig 5 and 6
writetable(coeff_results,'full_model_coeffs.xls','WriteRowNames',true)

%% results
load('results_zsubj.mat');

% behavioral results
orig_mean = nanmean(orig_distance); % in original paper units
orig_sd = nanstd(orig_distance);
irl_mean = nanmean(orig_distance)*6.55 % actual distance from target to dart "in real life"
irl_median = nanmedian(orig_distance)*6.55 % actual distance from target to dart "in real life"
irl_sd = nanstd(orig_distance)*6.55

orig_throwtime_mean = nanmean(orig_throwtime)
orig_throwtime_median = nanmedian(orig_throwtime)
orig_throwtime_sd = nanstd(orig_throwtime)

sum(r.pres)
sum(~r.pres)

% regression results
mdl_robust
sum(mdl_robust.Robust.Weights==0)
sortrows(mdl_robust.Residuals.Raw(mdl_robust.Robust.Weights==0))

% z = (sqrt(paper_dist)-mean(sqrt(paper_dist)))*std(sqrt(paper_dist))
z_errors = mdl_robust.Residuals.Raw; % in z units
t_unz_errors = (z_errors.*nanstd(z_distance)+nanmean(z_distance)); % transformed paper units
unz_errors = t_unz_errors.^2; % in paper units
irl_unz_errors = unz_errors*6.55; % at projection screen scale, cm, "in real life"
irl_mad_unz_error = nanmedian(abs(irl_unz_errors))

% figure 2a, exclude outliers for plotting sake
figure;
hist(mdl_robust.Residuals.Raw(mdl_robust.Robust.Weights~=0),100);
xlim([-5,5]);
saveas(gcf,'fig2a.jpg')

% figure 2b
figure;
y_hat = predict(mdl_robust, r);
plot(r.distance,y_hat,'bo'); axis([-5,5,-5,5]); hold on
plot([-5, 5], [-5,5])
saveas(gcf,'fig2b.jpg')

%% shuffle and split crossvalidation
nfolds = 1000;
training_mads = zeros(nfolds,1);
test_mads = zeros(nfolds,1);
training_rmses = zeros(nfolds,1);
test_rmses = zeros(nfolds,1);
adj_rs = zeros(nfolds,1);
training_resids = zeros(nfolds,1258);
training_inds = zeros(nfolds,1258);

test_resids = zeros(nfolds,139);
test_inds = zeros(nfolds,139);

for fold = 1:nfolds
	fold
	r.inds = (1:length(r.distance))';
	r.rand = randperm(length(r.distance))';
	r = sortrows(r, width(r));
	r = removevars(r,'rand');
	test_inds(fold,:)= r.inds(1:139);
	training_inds(fold,:)= r.inds(140:end);
	r = removevars(r,'inds');
	
	test_data = r(1:139,:);
	training_data = r(140:end,:);
	training_mdl = fitlm(training_data,'interactions','ResponseVar', 'distance', 'RobustOpts','on');
	training_resids(fold,:) = training_mdl.Residuals.Raw;
	test_resids(fold,:) = test_data.distance-predict(training_mdl, test_data);
	adj_rs(fold) = training_mdl.Rsquared.Adjusted;
end
training_mads = nanmedian(abs(training_resids),2);
test_mads = nanmedian(abs(test_resids),2);
mad_ratios = training_mads./test_mads;
save('cross_val_results.mat','adj_rs', 'training_resids','test_resids', 'training_mads',...
	'test_mads', 'mad_ratios');

%% permutation stats
nfolds = 1000;
rmses = nan(nfolds,1);
betas = nan(nfolds,379,4);
adj_rs = nan(nfolds,1);
resids = nan(nfolds,1397);
ps = nan(nfolds,1);
shuff_table = r;
for fold = 1:nfolds
	fold
	shuff_table.distance = shuff_table.distance(randperm(length(shuff_table.distance)));
	mdl = fitlm(shuff_table,'interactions','ResponseVar', 'distance', 'RobustOpts','on');
	res = table2array(anova(mdl,'summary'));
	ps(fold) = res(2,5);
	resids(fold,:) = mdl.Residuals.Raw;
	rmses(fold) = mdl.RMSE;
	betas(fold,:,:) = table2array(mdl.Coefficients);
	adj_rs(fold) = mdl.Rsquared.Adjusted;
end
mads = nanmedian(abs(resids),2);
save('rand_results.mat','adj_rs', 'resids','mads',...
	'rmses', 'betas','ps');

%% cross-validate figures
load cross_val_results
figure; % figure 3
hist(adj_rs,100);
mean(adj_rs)
median(adj_rs)
mean(adj_rs(adj_rs>.8))
saveas(gcf,'fig3.jpg')

figure; % figure 4
subplot(3,1,1); hist(training_mads,40); xlim([0 1.5])
subplot(3,1,2); hist(test_mads,40); xlim([0 1.5])
subplot(3,1,3); hist(mad_ratios,40); xlim([0 1.5])
saveas(gcf,'fig4.jpg')

sum(mdl_robust.Coefficients.pValue<0.05)
