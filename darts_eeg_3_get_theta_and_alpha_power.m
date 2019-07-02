clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/darts';
script_dir = 'C:\Users\Rob\Desktop\darts';
cd(script_dir);
rmpath('/data/common/matlab/eeglab')
addpath('./eeglab2019_0')
data_dir = './data/';
addpath(data_dir)
eeglab;
close all;

% user input
srate = 512;
chan_labs = {'Fz','Cz','Pz','Oz'};
freq_labs = {'theta', 'alpha', 'gamma'};
interval_labs = {'encoding','maintenance', 'recall'};
subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};

% load XLSX SNAP data
[num,txt,raw] = xlsread([data_dir,'behavioral_data.xlsx']);
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end

% initialize eeg columns to add to the behav columns
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs
			xlsx.([interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}]) = nan(length(xlsx.pres),1);
		end
	end
end

for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_',num2str(srate),'_laplace.set'];
	
	% bandpass filter at desired freqs
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	EEG_theta = pop_eegfiltnew(EEG, 3, 7);
	EEG_alpha = pop_eegfiltnew(EEG, 8, 12);
	EEG_gamma = pop_eegfiltnew(EEG, 21, 55);
	EEGs = {EEG_theta,EEG_alpha, EEG_gamma};
	load([data_dir, subj_id,'_eeg_',num2str(srate),'_latencys'])
	
	% correct data collection issues
	% the problem is there are some trials in the
	% xls file that are not in eeg
	n_snap_trials = sum(xlsx.subject == str2double(subj_id));
	n_eeg_trials = length(end_event_latencys);
	eeg_trial_strs = str2num(end_event_strings(:,1:4)); % ignore warning, use str2num
	subj_inds = xlsx.subject == str2double(subj_id);
	snap_trial_strs = str2num([num2str(xlsx.delay(subj_inds)),...
		num2str(xlsx.position(subj_inds),'%02d'),num2str(xlsx.pres(subj_inds))]);
	eeg_to_snap_inds = 1:length(eeg_trial_strs);
	if strcmp(subj_id, '580')
		eeg_to_snap_inds = 10+(1:n_eeg_trials);
	end
	if strcmp(subj_id,'616') || strcmp(subj_id,'621') || strcmp(subj_id,'627')
		eeg_to_snap_inds = [];
		for eeg_i = 1:length(eeg_trial_strs)
			for snap_i = eeg_i:length(snap_trial_strs)
				if eeg_trial_strs(eeg_i) == snap_trial_strs(snap_i)
					eeg_to_snap_inds = [eeg_to_snap_inds, snap_i];
					break
				end
			end
		end
	end
	eeg_to_snap_inds = eeg_to_snap_inds + find(xlsx.subject==str2double(subj_id),1) - 1;

	% retrieve this subjs trial-level amplitudes at desired chans and freqs
	for chan_lab_i = 1:length(chan_labs)
		chan_lab = chan_labs{chan_lab_i};
		for freq_lab_i = 1:length(EEGs)
			freq_lab = freq_labs{freq_lab_i};
			EEG = EEGs{freq_lab_i};
			EEG.data = abs(hilbert(EEG.data)).^2;
			maint_trial_amps = [];
			recall_trial_amps = [];
			
			for eeg_trial_i = 1:n_eeg_trials
				start_latency_i = start_event_latencys(eeg_trial_i);
				cue_latency_i = cue_event_latencys(eeg_trial_i);
				end_latency_i = end_event_latencys(eeg_trial_i);

				maint_baseline_inds = round(start_latency_i-0.2*srate):start_latency_i;
				maint_baseline = mean(EEG.data(chan_lab_i,maint_baseline_inds),2);
				maint_inds = start_latency_i:cue_latency_i;
				maint_trial = mean(EEG.data(chan_lab_i,maint_inds),2);
				maint_trial_amps(end+1) = maint_trial-maint_baseline;
				
				recall_baseline_inds = round(cue_latency_i-0.2*srate):cue_latency_i;
				recall_baseline = mean(EEG.data(chan_lab_i,recall_baseline_inds),2);
				recall_inds = cue_latency_i:end_latency_i;
				recall_trial = mean(EEG.data(chan_lab_i,recall_inds),2);
				recall_trial_amps(end+1) = recall_trial-recall_baseline;
			end
			
			maint_amps = xlsx.(['maintenance_',chan_lab,'_',freq_lab]);
			recall_amps = xlsx.(['recall_',chan_lab,'_',freq_lab]);
			for trial_amp_i = 1:n_eeg_trials
				maint_amps(eeg_to_snap_inds(trial_amp_i)) = maint_trial_amps(trial_amp_i);
				recall_amps(eeg_to_snap_inds(trial_amp_i)) = recall_trial_amps(trial_amp_i);
			end
			xlsx.(['maintenance_',chan_lab,'_',freq_lab]) = maint_amps;
			xlsx.(['recall_',chan_lab,'_',freq_lab]) = recall_amps;
	
		end
	end	
end

%% make regression vars
non_outlier_inds = (xlsx.distance < 6 & xlsx.throwtime < 5);
orig_throwtime = (xlsx.throwtime(non_outlier_inds));
orig_delay = xlsx.delaystilltime(non_outlier_inds);
pres = categorical(xlsx.pres(non_outlier_inds));
orig_dists = xlsx.distance(non_outlier_inds);
t_dists = orig_dists.^(1/2);

subjs = (xlsx.subject(non_outlier_inds));

% convert to z scores
delay = (orig_delay-nanmean(orig_delay))/nanstd(orig_delay);
dists = (t_dists-nanmean(t_dists))/nanstd(t_dists);
throwtime = (orig_throwtime-nanmean(orig_throwtime))/nanstd(orig_throwtime);

% create workspace variables for regression
chan_labs = {'Fz','Cz','Pz','Oz'};
freq_labs = {'theta', 'alpha','gamma'};
interval_labs = {'maintenance', 'recall'};
varnames = {'dists'; 'throwtime'; 'delay';'pres'}; % turn off for canonical corr
% varnames = {}; % turn on for canonical corr
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs
			varname = [interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}];
			varnames{end+1} = varname;
			eval([varname,' = (xlsx.',varname,'(non_outlier_inds)-nanmean(xlsx.',varname,'(non_outlier_inds)))/nanstd(xlsx.',varname,'(non_outlier_inds));'])
		end
	end
end
eval( ['reg_table = table(', strjoin(varnames,','),');'] );

%% behav+eeg regression
% eval( ['reg_table = table(', strjoin(varnames,','),');'] );
% mdl_ints = fitlm(reg_table,'interactions','ResponseVar', 'dists');
mdl_robust = fitlm(reg_table,'interactions','ResponseVar', 'dists', 'RobustOpts','on');
results = sortrows(anova(mdl_robust),5,'descend');
writetable(results,'full_model_results.xls','WriteRowNames',true)

%%
% plotResiduals(mdl_robust, 'Probability')
% plotResiduals(mdl_robust)
% mdl = stepwiselm(reg_table,'ResponseVar', 'dists')
% mdl = stepwiselm(reg_table,'interactions','ResponseVar', 'dists')
% 
% plot(reg_table.dists,predict(mdl_robust, reg_table),'.')
% 
% plot(reg_table.dists,mdl_robust.Residuals.Raw,'.')
% hist(mdl_robust.Residuals.Raw)
% results = sortrows(anova(mdl1),5,'descend');
nanmedian(abs(mdl_robust.Residuals.Raw))

%% leave-one-out crossvalidation
% loo_mads = [];
% loo_rmses = []; 
% loo_rs = [];
% eval( ['reg_table = table(', strjoin(varnames,','),');'] );
% for fold = 1:height(reg_table)
% 	fold
% 	eval( ['reg_table = table(', strjoin(varnames,','),');'] );
% 	reg_table(fold,:) = [];
% 	loo_mdl = fitlm(reg_table,'interactions','ResponseVar', 'dists', 'RobustOpts','on');
% 	loo_rs = [loo_rs; loo_mdl.Rsquared.Adjusted];
% 	loo_mads = [loo_mads; nanmedian(abs(loo_mdl.Residuals.Raw))];
% 	loo_rmses = [loo_rmses; loo_mdl.RMSE];
% end

%% shuffle and split crossvalidation
eval( ['reg_table = table(', strjoin(varnames,','),');'] );
nfolds = 10000;
training_mads = zeros(nfolds,1);
test_mads = zeros(nfolds,1);
training_rmses = zeros(nfolds,1);
test_rmses = zeros(nfolds,1);
adj_rs = zeros(nfolds,1);
training_resids = zeros(nfolds,1201);
training_inds = zeros(nfolds,1201);

test_resids = zeros(nfolds,133);
test_inds = zeros(nfolds,133);

for fold = 1:nfolds
	fold
	
	reg_table.inds = (1:length(reg_table.dists))';
	reg_table.rand = randperm(length(reg_table.dists))';
	reg_table = sortrows(reg_table, width(reg_table));
	reg_table = removevars(reg_table,'rand');
	test_inds(fold,:)= reg_table.inds(1:133);
	training_inds(fold,:)= reg_table.inds(134:end);
	reg_table = removevars(reg_table,'inds');

	test_data = reg_table(1:133,:);
	training_data = reg_table(134:end,:);
	training_mdl = fitlm(training_data,'interactions','ResponseVar', 'dists', 'RobustOpts','on');
	training_resids(fold,:) = training_mdl.Residuals.Raw;
	test_resids(fold,:) = test_data.dists-predict(training_mdl, test_data);
	adj_rs(fold) = training_mdl.Rsquared.Adjusted;

end
%%
high_test_inds = [];
high_training_inds = [];

for adj_r_i = 1:length(adj_rs)
	if adj_rs(adj_r_i) > 0.84
		high_training_inds = [high_training_inds, training_inds(adj_r_i,:)];
		high_test_inds = [high_test_inds, test_inds(adj_r_i,:)];

	end
end
%%
% histogram(high_training_inds,'BinMethod','integer')
[~, max_inds] = maxk(histcounts(high_training_inds,'BinMethod','integer'),134);
max_table = reg_table(max_inds,:);
max_table.pres = double(max_table.pres);
non_max_table = reg_table;
non_max_table(max_inds,:) = [];
tabulate(categorical(max_table.pres))
tabulate(categorical(max_table.subjs))
tabulate(categorical(reg_table.subjs))

ttest2(table2array(max_table),table2array(non_max_table))
%%
pres_inds = xlsx.pres(non_outlier_inds)==1;
y_hat = predict(mdl_robust, reg_table);
subplot(2,1,1);plot(reg_table.dists(~pres_inds),y_hat(~pres_inds),'ro'); axis(4*[-5,5,-5,5])
subplot(2,1,2);plot(reg_table.dists(pres_inds),y_hat(pres_inds),'bo'); axis(4*[-5,5,-5,5])
plot(reg_table.dists,y_hat,'bo'); axis([-5,5,-5,5]); hold on
plot([-5, 5], [-5,5])

irl_std = nanstd(orig_dists)*8.5;
orig_mean = nanmean(orig_dists);
orig_sd = nanstd(orig_dists);
z_errors = mdl_robust.Residuals.Raw; % in z units
unz_errors = ((z_errors+orig_mean).*orig_sd).^2; % in original paper units
irl_unz_errors = unz_errors*8.5; % irl units

median_z_error = nanmedian(abs(z_errors));
median_unz_error = nanmedian(abs(unz_errors));
irl_median_unz_error = median_unz_error*8.5; % convert paper units to irl units

median(nanmedian(abs(training_resids),2))
hist(nanmedian(abs(training_resids),2))
hist(nanmedian(abs(test_resids),2))

training_mads = nanmedian(abs(training_resids),2);
test_mads = nanmedian(abs(test_resids),2);
mad_ratios = training_mads./test_mads;
figure;
subplot(3,1,1); hist(training_mads,30); xlim([0 1.5])
subplot(3,1,2); hist(test_mads,30); xlim([0 1.5])
subplot(3,1,3); hist(mad_ratios,30); xlim([0 1.5])

	% median absolute difference
% 	training_mads(fold) = nanmedian(abs(training_mdl.Residuals.Raw));
% 	test_mads(fold) = nanmedian(abs(test_resids(fold,:)));
	
	% root-mean-square error
% 	training_rmses(fold) = training_mdl.RMSE;
% 	test_rmses(fold) = sqrt(nanmean(test_resids(fold,:).^2));

%% canonical correlation analysis
% X = [dists,throwtime,delay, pres];
% for col = 1:size(X,2)
% 	col_mean = nanmean(X(:,col));
% 	for row = 1:size(X,1)
% 		if isnan(X(row,col))
% 			X(row,col) = col_mean;
% 		end
% 	end
% end
% eval( ['Y = [', strjoin(varnames,','),'];'] );
% for col = 1:size(Y,2)
% 	col_mean = nanmean(Y(:,col));
% 	for row = 1:size(Y,1)
% 		if isnan(Y(row,col))
% 			Y(row,col) = col_mean;
% 		end
% 	end
% end
% 
%  [A,B,R,U,V,STATS] = canoncorr(X,Y);
%  figure
%  plot(U(:,1),V(:,1),'.')
%  [rho,pval] = corr(U(:,1),V(:,1), 'type','Pearson');
%  canon_variate_x = A(:,1);
%  canon_variate_y = B(:,1);
% 
%  plotnames_x = categorical({'dists'; 'throwtime'; 'delay';'pres'});
% 
%  bar(plotnames_x,canon_variate_x)
%  plotnames = {};
%  for chan_lab = chan_labs
% 	for freq_lab = freq_labs
% 		for interval_lab = interval_labs
% 			plotname = [interval_lab{:},' ',chan_lab{:},' ',freq_lab{:}];
% 			plotnames{end+1} = plotname;
% 		end
% 	end
%  end
%  [sorted_can_var_y, idx] = sort(canon_variate_y, 'descend','ComparisonMethod','abs');
%  for i = 1:length(idx)
% 	 sorted_plotnames(i) = plotnames(idx(i));
%  end
%   c = categorical(sorted_plotnames);
% 	c = reordercats(c,sorted_plotnames);
% 
%  bar(c, sorted_can_var_y)
%   [rho,pval]

%%
% fcn = @(reg_table) predict( fitlm(reg_table,'interactions','ResponseVar', 'dists','RobustOpts','on'), reg_table.dists);
% perform cross-validation, and return average MSE across folds
% mse = crossval('mse', X, Y,'Predfun',fcn, 'kfold',10);
% compute root mean squared error
% avrg_rmse = sqrt(mse)
% results_array = table2array(results);
% sigs = results.Properties.RowNames(find(results_array(:,5)<0.05));
% mdl2 = fitlm(reg_table, ['dists ~ ',strjoin(sigs,' + ')],'RobustOpts','on')

% plotResiduals(fit_eeg, 'Probability')

% mdl = fitlm(reg_table,'interactions','ResponseVar', 'dists', 'RobustOpts', 'on')
% plotResiduals(mdl, 'Probability')
% 
% fit_eeg = stepwiselm(reg_table,'interactions','ResponseVar', 'dists')
% % larg = find(abs(mdl.Residuals.Raw)>4);
% larg = [998, 1282, 1321, 1285, 1313]
% mdl2 = fitlm(reg_table,'interactions','ResponseVar', 'dists', 'RobustOpts', 'on','Exclude',larg)
% plotResiduals(mdl2, 'Probability')
% 
% mdl3 = fitlm(reg_table,'interactions','ResponseVar', 'dists','Exclude',larg)
% 
% [~,larg] = max(mdl.Diagnostics.CooksDistance,10);
% plot(abs(mdl.Residuals.Raw),'o')
% plotResiduals(mdl)
% plotResiduals(mdl, 'Probability')
% 
% plotDiagnostics(mdl)
% figure; plotDiagnostics(mdl,'cookd')
% plotResiduals(mdl,'lagged')
% plotResiduals(mdl,'fitted')
% step(mdl)

%% correlation
% eval( ['array = [', strjoin(varnames,','),'];'] );
% nans = sum(isnan(array),2) > 0;
% array = array(~nans,:);
% [rho,pval] = corr(array, 'rows','complete');
% heatmap(rho)
% heatmap(pval)

%% partial least squares regression
% eval( ['array = [', strjoin(varnames,','),'];'] );
% nans = sum(isnan(array),2) > 0;
% array = array(~nans,:);
% [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]= plsregress(array(:,2:end),array(:,1),size(array,2)-1);
% plot(1:length(PCTVAR),cumsum(100*PCTVAR(2,:)),'-bo');

%% test collinearity
% col_varnames = varnames;
% col_varnames(ismember(col_varnames,{'dists','pres'})) = [];
% col_reg_table = removevars(reg_table,{'dists','pres'});
% collintest(col_reg_table)
% corrplot(col_reg_table,'testR','on')
% VIF = diag(inv(corr(table2array(col_reg_table))))'
% results = sortrows(anova(fit_eeg),5,'descend')
% results = table2array(anova(fit_eeg));
% results(results(:,5) < 0.05)
% mdl = stepwiselm(reg_table,'interactions','ResponseVar', 'dists')
% 
% thresholds = 0.05./(height(results):-1:1)'
% sum(results.pValue<thresholds)
% fdr_res = fdr_bh(results.pValue, 0.9999, 'dep', 'yes');
% for k = 1:length(results.pValue)
% 	if results.pValue(k) <= (k/length(results.pValue))*0.05
% 		break
% 	end
% end
% m = length(results.pValue);
% results.pValue<=0.05
% plot(results.pValue)
%%
% for varname = varnames(4:end)
% 	varnames{end+1} = [varname{1};'*throwtime']
% 	varnames{end+1} = [varname{1},'*delay']
% end
% 
% mdl = fitlm(reg_table, ['dists ~ ',strjoin(varnames(2:end),' + ')], 'RobustOpts', 'on')
% mdl = stepwiselm(reg_table, ['dists ~ ',strjoin(varnames(2:end),' + ')])

%% behavioral regression
% varnames = {'dists', 'throwtime', 'delay', 'pres'};
% eval( ['reg_table = table(', strjoin(varnames,','),');'] );
% fit_behav = fitlm(reg_table,'interactions', 'ResponseVar', 'dists', 'RobustOpts', 'on')
% swfit_behav = stepwiselm(reg_table,'interactions', 'ResponseVar', 'dists', 'Criterion', 'AdjRsquared')
% anova(swfit_behav)