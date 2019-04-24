% Script timing
clear
clc
start = datetime('now');
debugging = 1;

% Initialize paths
warning('off','MATLAB:rmpath:DirNotFound')
warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:MKDIR:DirectoryExists')
scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
cd(scriptdir)
rmpath /data/common/matlab/eeglab/
addpath(genpath('./eeglab14_1_2b'))

if debugging == 1
	inpath = './coupled_osc_32/';
	cd(inpath)
	mats = dir('*.mat');
	cd(scriptdir)
	outpath = './corr_32/';
	mkdir(outpath)
	srate = 32;
else
	inpath = './coupled_osc_512/';
	cd(inpath)
	mats = dir('*.mat');
	cd(scriptdir)
	outpath = './corr_512/';
	mkdir(outpath)
	srate = 512;
end

% load XLSX SNAP data
[num,txt,raw] = xlsread('/data/mobi/Darts/Data/Analysis Data/XLSX/Combined.xlsx');
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end
subjs = unique(xlsx.subject)';

% for each subj_Hz
for subj_i = 1:length(subjs)
	subj_n = subjs(subj_i);
	subj_id = num2str(subj_n);
	subj_ps = [];
	subj_rs = [];
	prctdone(subj_i-1,length(subjs));
	
	% load trial indices
	trialpath = './trial_inds_32/';
	trialmats = dir([trialpath,'*.mat']);
	for trialmat_i = 1:length(trialmats)
		if strcmp(trialmats(trialmat_i).name(18:20),subj_id)
			try
				trial_inds_mat = load([trialpath,trialmats(trialmat_i).name]);
			catch
				error('Could not find trial timings')
			end
		end
	end
	eeg_trial_strs = str2num(trial_inds_mat.end_event_strings(:,1:4)); %#ok<*ST2NM>
	end_event_time_inds = trial_inds_mat.end_event_time_inds;
	start_event_time_inds = trial_inds_mat.start_event_time_inds;
	n_time_points = max(end_event_time_inds-start_event_time_inds)-srate/2;
	
	% correct eeg to snap index mapping
	n_eeg_trials = length(end_event_time_inds);
	subj_inds = xlsx.subject == subj_n;
	snap_trial_strs = str2num([num2str(xlsx.delay(subj_inds)),num2str(xlsx.position(subj_inds),'%02d'),num2str(xlsx.pres(subj_inds))]);
	n_snap_trials = length(snap_trial_strs);
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
	
	% load TO and FROM spec data
	for to_Hz = 1:15
		
		% load TO spec data
		specpath = './hilbert_ics_32/';
		specmats = dir([specpath,'*.mat']);
		to_spec_data = [];
		for to_specmat_i = 1:length(specmats)
			if strcmp(specmats(to_specmat_i).name(18:20),subj_id)
				if str2double(specmats(to_specmat_i).name(29:31)) == to_Hz
					to_spec_data_mat = load([specpath,specmats(to_specmat_i).name]);
					to_spec_data = double(to_spec_data_mat.spec_data);
					break
				end
			end
		end

		% load FROM coupling and spec data
		for mat_i = 1:length(mats)
			if ~strcmp(subj_id,mats(mat_i).name(1:3)) || ~strcmp(num2str(to_Hz), mats(mat_i).name(5:(end-4)))
				continue
			end
			load([inpath,mats(mat_i).name]);
			from_Hzs = [from_Hzs(1), from_Hzs(2), from_Hzs(2), from_Hzs(3)];
			from_props = [from_props(1), from_props(2), from_props(2), from_props(3)];
			from_spec_data = [];
			for to_prop_i = 1:4
				from_Hz = from_Hzs(to_prop_i);
				from_prop = from_props(to_prop_i);
				spec_prop_temp = 'APPF';
				from_prop_i = spec_prop_temp == from_prop;
				for from_specmat_i = 1:length(specmats)
					if strcmp(specmats(from_specmat_i).name(18:20),subj_id)
						if str2double(specmats(from_specmat_i).name(29:31)) == from_Hz
							from_spec_data_mat = load([specpath,specmats(from_specmat_i).name]);
							from_spec_data(to_prop_i,:) = double(from_spec_data_mat.spec_data(from_prop_i,:));
							break
						end
					end
				end
			end
		end
		
		% get trial behav and spec data
		n_props = size(to_spec_data,1);
		subj_data = structfun(@(F) F([xlsx.subject]==subj_n), xlsx, 'uniform', 0);
		trial_dists = [];
		trial_sqrt_dists = [];
		trial_2_dists = [];
		trial_preses = [];
		trial_amps = [];
		trial_freqs = [];
		trial_resids = [];
		warning('off','MATLAB:rankDeficientMatrix')
		for eeg_trial_i = 1:n_eeg_trials
			% get behav data
			if eeg_trial_strs(eeg_trial_i)~=snap_trial_strs(eeg_to_snap_inds(eeg_trial_i))
				error('Trial mismatch')
			end
			trial_dist = subj_data.distance(eeg_to_snap_inds(eeg_trial_i));
			trial_sqrt_dist = trial_dist.^(1/2);
			trial_2_dist = trial_dist.^2;
			trial_pres = subj_data.pres(eeg_to_snap_inds(eeg_trial_i));
			if trial_2_dist> 30
				continue
			end
			trial_dists = [trial_dists; trial_dist];
			trial_sqrt_dists = [trial_sqrt_dists; trial_sqrt_dist];
			trial_2_dists = [trial_2_dists; trial_2_dist];
			trial_preses = [trial_preses; trial_pres];
			
			% get spec data
			to_trial_data = to_spec_data(1:n_props,start_event_time_inds(eeg_trial_i):start_event_time_inds(eeg_trial_i)+n_time_points)';
			trial_Y = detrend(to_trial_data-mean(to_trial_data)); % detrends and centers
			from_trial_data = from_spec_data(1:n_props,start_event_time_inds(eeg_trial_i):start_event_time_inds(eeg_trial_i)+n_time_points)';
			trial_X = detrend(from_trial_data-mean(from_trial_data)); % detrends and centers
			trial_amps = [trial_amps;mean(trial_Y(1,:))];
			trial_freqs = [trial_freqs;mean(trial_Y(4,:))];
			
			% CFC linear model, Y = XB => X^(-1)*Y = B
			trial_B = trial_X\trial_Y;
			trial_resid = min(mean((trial_Y-trial_X*trial_B).^2));
			trial_resids = [trial_resids; trial_resid];
		end
		
		% p-values
		is = logical(trial_preses);
		[r,p] = corr(trial_dists(is),trial_amps(is),'rows','complete','type','pearson');
		subj_pres_amp_ps(to_Hz,:) = p;
		[r,p] = corr(trial_dists(is),trial_freqs(is),'rows','complete','type','pearson');
		subj_pres_freq_ps(to_Hz,:) = p;
		[r,p] = corr(trial_dists(is),trial_resids(is,:),'rows','complete','type','pearson');
		subj_pres_resid_ps(to_Hz,:) = p;

		is = ~logical(trial_preses);
		[r,p] = corr(trial_dists(is),trial_amps(is),'rows','complete','type','pearson');
		subj_abs_amp_ps(to_Hz,:) = p;
		[r,p] = corr(trial_dists(is),trial_freqs(is),'rows','complete','type','pearson');
		subj_abs_freq_ps(to_Hz,:) = p;
		[r,p] = corr(trial_dists(is),trial_resids(is,:),'rows','complete','type','pearson');
		subj_abs_resid_ps(to_Hz,:) = p;
		
		% r values
		is = logical(trial_preses);
		[r,p] = corr(trial_dists(is),trial_amps(is),'rows','complete','type','pearson');
		subj_pres_amp_rs(to_Hz,:) = r;
		[r,p] = corr(trial_dists(is),trial_freqs(is),'rows','complete','type','pearson');
		subj_pres_freq_rs(to_Hz,:) = r;
		[r,p] = corr(trial_dists(is),trial_resids(is,:),'rows','complete','type','pearson');
		subj_pres_resid_rs(to_Hz,:) = r;

		is = ~logical(trial_preses);
		[r,p] = corr(trial_dists(is),trial_amps(is),'rows','complete','type','pearson');
		subj_abs_amp_rs(to_Hz,:) = r;
		[r,p] = corr(trial_dists(is),trial_freqs(is),'rows','complete','type','pearson');
		subj_abs_freq_rs(to_Hz,:) = r;
		[r,p] = corr(trial_dists(is),trial_resids(is,:),'rows','complete','type','pearson');
		subj_abs_resid_rs(to_Hz,:) = r;
	end
	% p-values
	all_pres_amp_ps(:,subj_i) = subj_pres_amp_ps;
	all_pres_freq_ps(:,subj_i) = subj_pres_freq_ps;
 	all_pres_resid_ps(:,subj_i) = subj_pres_resid_ps;
	all_abs_amp_ps(:,subj_i) = subj_abs_amp_ps;
	all_abs_freq_ps(:,subj_i) = subj_abs_freq_ps;
 	all_abs_resid_ps(:,subj_i) = subj_abs_resid_ps;
	% r-values
	all_pres_amp_rs(:,subj_i) = subj_pres_amp_rs;
	all_pres_freq_rs(:,subj_i) = subj_pres_freq_rs;
 	all_pres_resid_rs(:,subj_i) = subj_pres_resid_rs;
	all_abs_amp_rs(:,subj_i) = subj_abs_amp_rs;
	all_abs_freq_rs(:,subj_i) = subj_abs_freq_rs;
 	all_abs_resid_rs(:,subj_i) = subj_abs_resid_rs;	
end

	all_abs_zs = .5*(log(1+all_abs_resid_rs)-log(1-all_abs_resid_rs));
	all_pres_zs = .5*(log(1+all_pres_resid_rs)-log(1-all_pres_resid_rs));
	
[H,P,CI,STATS] = ttest(all_abs_zs', 0);
all_abs_mean = mean(all_abs_zs');
all_abs_sem = std(all_abs_zs')/sqrt(length(all_abs_zs'));

%% Script timing
stop = datetime('now');
disp(['Start ',datestr(start)])
disp(['End ',datestr(stop)])
disp(['Runtime ', char(duration(stop-start))])