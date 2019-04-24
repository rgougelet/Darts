%% Script timing
clear
start = datetime('now');
debugging = 1;

%% Initialize paths
warning('off','MATLAB:rmpath:DirNotFound')
scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
cd(scriptdir)
rmpath /data/common/matlab/eeglab/
addpath(genpath('./eeglab14_1_2b'))

if debugging == 1
	inpath = './hilbert_ics_32/';
	cd(inpath)
	mats = dir('*.mat');
	cd(scriptdir)
	outpath = './reg_mats_32/';
	mkdir(outpath)
	srate = 32;
else
	inpath = './hilbert_ics_512/';
	cd(inpath)
	mats = dir('*.mat');
	cd(scriptdir)
	outpath = './reg_mats_512/';
	mkdir(outpath)
	srate = 512;
end

for subj_id = {'580','627','616','608','579','621','619','607','631','571'}
	subj_set_inds = [];
	for comp_freq_ind = 1:length(mats)
		if strcmp(mats(comp_freq_ind).name(18:20),subj_id{:})
			subj_set_inds=[subj_set_inds,comp_freq_ind];
		end
	end
	
	%% Create "oscillation" regression matrix for each of this subject's "oscillations"
	reg_mat = [];
	for comp_freq_ind = subj_set_inds
		
		% load spectral data
% 		disp(['Appending regression matrix with ',mats(comp_freq_ind).name]);
		spec_data_mat = load([inpath,mats(comp_freq_ind).name]);
		spec_data = double(spec_data_mat.spec_data);
		comp_entropy = spec_data_mat.comp_entropy;
		n_props = size(spec_data,1);
		
		% import trial start and end indices
		try
			trial_time_mat = load([scriptdir,'/trial_inds_32/',mats(comp_freq_ind).name(1:20),'_AD_32.mat']);
		catch
			error('Could not find trial timings')
		end
		end_event_time_inds = trial_time_mat.end_event_time_inds;
		start_event_time_inds = trial_time_mat.start_event_time_inds;
		n_time_points = max(end_event_time_inds-start_event_time_inds)-srate/2;
		n_trials = length(end_event_time_inds);
		
		% localize to only within trials.
		osc_spec_data = nan(n_props,n_time_points,n_trials);
		for trial_i = 1:n_trials
					trial_data = spec_data(:,start_event_time_inds(trial_i):end_event_time_inds(trial_i)-srate/2)';
% 			trial_data = spec_data(1:n_props,start_event_time_inds(trial_ind):start_event_time_inds(trial_ind)+n_time_points)';
			detrended_trial_data = detrend(trial_data-mean(trial_data)); % detrends and centers
			osc_spec_data(1:n_props,1:length(trial_data),trial_i) = detrended_trial_data';
		end
		reg_mat = cat(4,reg_mat, osc_spec_data);
	end
	disp(['Regression matrix ',num2str(100*sum(isnan(reg_mat(:)))/numel(reg_mat)),' percent NaN'])
	parsave([outpath, subj_id{:},'.mat'],{reg_mat},{'reg_mat'});
end

%% Script timing
stop = datetime('now');
disp(['Start ',datestr(start)])
disp(['End ',datestr(stop)])
disp(['Runtime ', char(duration(stop-start))])