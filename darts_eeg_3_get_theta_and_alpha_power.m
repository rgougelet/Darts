clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/darts';
script_dir = 'C:\Users\Rob\Desktop\darts';
cd(script_dir);
rmpath('/data/common/matlab/eeglab')
addpath('./eeglab2019_0')
data_dir = './data/';
addpath(data_dir)

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
	EEG_gamma = pop_eegfiltnew(EEG, 21, 31);
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
% 			encode_trial_amps = [];
			maint_trial_amps = [];
			recall_trial_amps = [];
			for eeg_trial_i = 1:n_eeg_trials
				start_latency_i = start_event_latencys(eeg_trial_i);
				cue_latency_i = cue_event_latencys(eeg_trial_i);
				end_latency_i = end_event_latencys(eeg_trial_i);
% 				encode_trial_amps(end+1) = mean(abs(EEG.data(chan_lab_i,start_latency_i:(start_latency_i+srate))),2)';

				maint_baseline_inds = round(start_latency_i-0.2*srate):start_latency_i;
% 				maint_baseline = sum(EEG.data(chan_lab_i,maint_baseline_inds),2);
				maint_baseline = mean(EEG.data(chan_lab_i,maint_baseline_inds),2);
% 				maint_baseline = mean(EEG.data(chan_lab_i,maint_baseline_inds),2)/length(maint_baseline_inds);

				maint_inds = start_latency_i:cue_latency_i;
% 				maint_trial = sum(EEG.data(chan_lab_i,maint_inds),2);
				maint_trial = mean(EEG.data(chan_lab_i,maint_inds),2);
% 				maint_trial = sum(EEG.data(chan_lab_i,maint_inds),2)/length(maint_inds);
				
				maint_trial_amps(end+1) = maint_trial-maint_baseline;
				
				recall_baseline_inds = round(cue_latency_i-0.2*srate):cue_latency_i;
% 				recall_baseline = sum(EEG.data(chan_lab_i,recall_baseline_inds),2);
				recall_baseline = mean(EEG.data(chan_lab_i,recall_baseline_inds),2);
% 				recall_baseline = sum(EEG.data(chan_lab_i,recall_baseline_inds),2)/length(recall_baseline_inds);
				
				recall_inds = cue_latency_i:end_latency_i;
% 				recall_trial = sum(EEG.data(chan_lab_i,recall_inds),2);
				recall_trial = mean(EEG.data(chan_lab_i,recall_inds),2);
% 				recall_trial = sum(EEG.data(chan_lab_i,recall_inds),2)/length(recall_inds);
				
				recall_trial_amps(end+1) = recall_trial-recall_baseline;

			end
			
% 			encode_chan_freq_amps = xlsx.(['encoding_',chan_lab,'_',freq_lab]);
			maint_amps = xlsx.(['maintenance_',chan_lab,'_',freq_lab]);
			recall_amps = xlsx.(['recall_',chan_lab,'_',freq_lab]);

			for trial_amp_i = 1:n_eeg_trials
% 				encode_chan_freq_amps(eeg_to_snap_inds(trial_amp_i)) = encode_trial_amps(trial_amp_i);
				maint_amps(eeg_to_snap_inds(trial_amp_i)) = maint_trial_amps(trial_amp_i);
				recall_amps(eeg_to_snap_inds(trial_amp_i)) = recall_trial_amps(trial_amp_i);
			end
% 			xlsx.(['encoding_',chan_lab,'_',freq_lab]) = encode_chan_freq_amps;
% 			maint_amps(isnan(maint_amps))=nanmean(maint_amps);
% 			recall_amps(isnan(recall_amps))=nanmean(recall_amps);
			xlsx.(['maintenance_',chan_lab,'_',freq_lab]) = maint_amps;
			xlsx.(['recall_',chan_lab,'_',freq_lab]) = recall_amps;
	
		end
	end	
end

%% make regression vars
non_outlier_inds = (xlsx.distance < 6 & xlsx.throwtime < 5);
throwtime = xlsx.throwtime(non_outlier_inds);
delay = xlsx.delaystilltime(non_outlier_inds);
pres = categorical(xlsx.pres(non_outlier_inds));
dists = xlsx.distance(non_outlier_inds);
subjs = categorical(xlsx.subject(non_outlier_inds));
delay = (delay-nanmean(delay))/nanstd(delay);
dists = (dists-nanmean(dists))/nanstd(dists);
throwtime = (throwtime-nanmean(throwtime))/nanstd(throwtime);

chan_labs = {'Fz','Cz','Pz','Oz'};
freq_labs = {'theta', 'alpha','gamma'};
interval_labs = {'maintenance', 'recall'};
varnames = {'dists'; 'throwtime'; 'delay';'pres'};
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs
			varname = [interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}];
			varnames{end+1} = varname;
% 			eval([varname,' = xlsx.',varname,'(non_outlier_inds);'])
			eval([varname,' = (xlsx.',varname,'(non_outlier_inds)-nanmean(xlsx.',varname,'(non_outlier_inds)))/nanstd(xlsx.',varname,'(non_outlier_inds));'])
% 			eval(['figure;hist(',varname,')'])
		end
	end
end

% behav+eeg regression
eval( ['reg_table = table(', strjoin(varnames,','),');'] );
% dists = 1+(1:2000)';
% dists2 = 3*(2:2001)';
% reg_table = table(dists,dists2)
mdl_basic = fitlm(reg_table,'ResponseVar', 'dists');
mdl_ints = fitlm(reg_table,'interactions','ResponseVar', 'dists');

mdl_robust = fitlm(reg_table,'interactions','ResponseVar', 'dists', 'RobustOpts','on');
% mdl1 = fitlm(reg_table,'dists ~ recall_Cz_gamma:recall_Pz_gamma + maintenance_Oz_theta:maintenance_Oz_gamma', 'RobustOpts','on')
% mdl = stepwiselm(reg_table,'interactions','ResponseVar', 'dists')
% sw_mdl_r_squared = stepwiselm(reg_table,'interactions','ResponseVar', 'dists', 'Criterion', 'AdjRsquared')

% plotResiduals(mdl1, 'Probability')
% results = sortrows(anova(mdl1),5,'descend')

% hist(mdl1.Residuals.Raw)

%%
mad_ratios = [];
rmse_ratios = [];
rs = [];
for fold = 1:height(reg_table)
	fold
	eval( ['reg_table = table(', strjoin(varnames,','),');'] );
% 	reg_table.rand = randperm(length(reg_table.dists))';
% 	reg_table = sortrows(reg_table, width(reg_table));
% 	reg_table = removevars(reg_table,'rand');
	reg_table(fold,:) = [];
% 	test_data = reg_table(1:133,:);
% 	training_data = reg_table(134:end,:);
	loo_mdl = fitlm(reg_table,'interactions','ResponseVar', 'dists');
	rs = [rs; loo_mdl.Rsquared.Adjusted];
% 	hist(reg_table.dists)
% 	hist(training_mdl.Residuals.Raw)
% 		hist(test_resids)

	training_mad = nanmedian(abs(training_mdl.Residuals.Raw));
	test_resids = test_data.dists-predict(training_mdl, test_data);
	test_mad = nanmedian(abs(test_resids));
	mad_ratios = [mad_ratios; training_mad/test_mad];
	
	training_rmse = training_mdl.RMSE;
	test_rmse = sqrt(nansum(test_resids).^2);
	rmse_ratios = [rmse_ratios; training_rmse/test_rmse];
% 	y = nansum((test_data.dists-nanmean(test_data.dists)).^2);
% 	g = nansum((test_data.dists-predict(mdl1, test_data)).^2);
% 	r_squared = 1 - (g/y)
% 	mdl1.mse
% 	rs = [rs; corr(test_data.dists,predict(mdl1, test_data), 'rows', 'complete', 'type', 'Pearson')];

% 	rhos = [rhos; corr(test_data.dists,predict(mdl1, test_data), 'rows', 'complete', 'type', 'Spearman')];
% 	sqrt(nanmean((test_data.dists-predict(mdl1, test_data)).^2))
% 	nanstd(test_data.dists)
end
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

%% canonical correlation
% X = array(:,1:4);
% Y = array(:,5:end);
% [A,B,r,U,V] = canoncorr(X,Y)
% plot(U(:,1),V(:,1),'.')

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


% % [H,P,CI,STATS] = ttest2(Cz_thetas_abs,Cz_thetas_pres)
% % [H,P,CI,STATS] = ttest2(Cz_thetas_abs,Cz_thetas_pres, 'vartype', 'unequal')

%% behavioral regression
varnames = {'dists', 'throwtime', 'delay', 'pres'};
eval( ['reg_table = table(', strjoin(varnames,','),');'] );
fit_behav = fitlm(reg_table,'interactions', 'ResponseVar', 'dists', 'RobustOpts', 'on')
% swfit_behav = stepwiselm(reg_table,'interactions', 'ResponseVar', 'dists', 'Criterion', 'AdjRsquared')
% anova(swfit_behav)