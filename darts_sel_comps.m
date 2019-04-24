%% Script timing
clear
start = datetime('now');
debugging = 1;

%% Initialize paths
warning('off','MATLAB:rmpath:DirNotFound')
scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
cd(scriptdir)
rmpath /data/common/matlab/eeglab/
rmpath(genpath('./eeglab14_1_2b'))
addpath('./eeglab14_1_2b')
eeglab
if debugging == 1
	inpath = './eeg_32/';
	cd(inpath)
	sets = dir('*.set');
	cd(scriptdir)
	outpath = './aac_32/';
	mkdir(outpath)
	srate = 32;
else
	inpath = './eeg_512/';
	cd(inpath)
	sets = dir('*.set');
	cd(scriptdir)
	outpath = './aac_512/';
	mkdir(outpath)
end

%% For each subject:freq
for set_i = 1:length(sets)
	prctdone(set_i-1,length(sets));
	
	% load data
	eeg1= pop_loadset('filename', sets(set_i).name, 'filepath', inpath);
	eeg1 = eeg_checkset( eeg1 );
	eeg1.icaact = (eeg1.icaweights*eeg1.icasphere)*eeg1.data(eeg1.icachansind,:);
	theta_EEG = pop_eegfiltnew(eeg1, 3,7,[],0,[],0);
	alpha_EEG = pop_eegfiltnew(eeg1, 8,12,[],0,[],0);
	theta_EEG.icaact = (theta_EEG.icaweights*theta_EEG.icasphere)*theta_EEG.data(theta_EEG.icachansind,:);
	alpha_EEG.icaact = (alpha_EEG.icaweights*alpha_EEG.icasphere)*alpha_EEG.data(alpha_EEG.icachansind,:);
	subj_id = eeg1.setname(18:20);
	
	% import trial start and end indices
	try
		trial_time_mat = load([scriptdir,'/trial_inds_32/',eeg1.setname(1:20),'_AD_32.mat']);
	catch
		error('Could not find trial timings')
	end
	end_event_time_inds = trial_time_mat.end_event_time_inds;
	start_event_time_inds = trial_time_mat.start_event_time_inds;
	n_time_points = max(end_event_time_inds-start_event_time_inds)-srate/2;
	n_trials = length(end_event_time_inds);
	
	% trim to trial data
	theta_trial_data = [];
	alpha_trial_data = [];
	for trial_i = 1:length(start_event_time_inds)
		start_i = start_event_time_inds(trial_i);
		end_i = end_event_time_inds(trial_i);
		theta_trial_data = cat(2,theta_trial_data,theta_EEG.icaact(:,start_i:end_i));
		alpha_trial_data = cat(2,alpha_trial_data,alpha_EEG.icaact(:,start_i:end_i));
	end

	% component x component aac_glm
	% hi = abs(hilbert(alpha_trial_data'));
	% lo = abs(hilbert(theta_trial_data'));

	% aacs = NaN(length(eeg1.icachansind),length(eeg1.icachansind));
	% for c1 = 1:length(eeg1.icachansind)
		% prctdone(c1,length(eeg1.icachansind)+1);
		% for c2 = 1:length(eeg1.icachansind)
			% y = (hi(:,c2)-mean(hi(:,c2)))/std(hi(:,c2));
			% X = (lo(:,c1)-mean(lo(:,c1)))/std(lo(:,c1));
			% [b, bint, resid] = regress(y, X);
			% aacs(c1,c2) = 1 - (resid'*resid)/(y'*y);
		% end
	% end
	% save output
	% parsave([sets(set_i).name(18:20), '_aac'],aacs,'aacs')

	% component x component pac_glm
	hi = abs(hilbert(alpha_trial_data'));
	lo = angle(hilbert(theta_trial_data'));
	Xc = cos(lo);
	Xs = sin(lo);

	pacs = NaN(length(eeg1.icachansind),length(eeg1.icachansind));
	for c1 = 1:length(eeg1.icachansind)
		prctdone(c1,length(eeg1.icachansind)+1);
		for c2 = 1:length(eeg1.icachansind)
			y = (hi(:,c2)-mean(hi(:,c2)))/std(hi(:,c2));
			X = [(Xc(:,c1)-mean(Xc(:,c1)))/std(Xc(:,c1)), (Xs(:,c1)-mean(Xs(:,c1)))/std(Xs(:,c1))];
			[b, bint, resid] = regress(y, X);
			pacs(c1,c2) = 1 - (resid'*resid)/(y'*y);
		end
	end
	% save output
	parsave([sets(set_i).name(18:20), '_pac'],pacs,'pacs')

end
