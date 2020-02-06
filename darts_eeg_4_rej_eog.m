clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
data_dir = [script_dir,'data/'];
addpath(data_dir)
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
eeglab;

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
%%
for subj_i = 1:length(subjs_to_include)
	close all;

	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'*_ic.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	% get eog inds
	uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
	lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
	lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
	rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
	veog = EEG.data(uveog_i,:)-EEG.data(lveog_i,:);
	heog = EEG.data(lheog_i,:)-EEG.data(rheog_i,:);

	% identify rej components
	veog_rej_inds = abs(corr(veog',EEG.icaact(:,:)'))>0.08;
	heog_rej_inds = abs(corr(heog',EEG.icaact(:,:)'))>0.08;
	eog_rej_inds = find(heog_rej_inds | veog_rej_inds);

	% rej eog ics
	EEG = pop_subcomp( EEG, eog_rej_inds, 0);
	EEG.etc.pipeline{end+1} =  ['ICs rejected: ', num2str(eog_rej_inds)];
	EEG = eeg_checkset(EEG, 'ica');

	[icx,f] = pwelch(EEG.icaact(:,:)',EEG.srate*100,EEG.srate*50,[],EEG.srate,'onesided');
	plot(f,icx)
	continue

	[~, ltheta_i] = min(abs(3-f));
	[~, rtheta_i] = min(abs(8-f));
	theta_i = [ltheta_i, rtheta_i];
	[~, ralpha_i] = min(abs(12-f));
	alpha_i = [rtheta_i+1, ralpha_i];
	[~, llogamma_i] = min(abs(30-f));
	[~, rlogamma_i] = min(abs(8-f));
	theta_i = [ltheta_i, rtheta_i];


	% remove eog
	EEG = pop_select( EEG,'nochannel', {'UVEOG', 'LVEOG', 'LHEOG', 'RHEOG'});
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'EOG channels removed';

	% align head model, warp to fiducials and Cz(45)->C21(85), Iz(88)->D23(119)
	EEG = headfit(EEG,subj_id);

	ic_EEG = EEG;

	% load pre-ica data
	subj_set = dir([data_dir,subj_id,'_eeg_512.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;

	% reject channels
	if strcmp(subj_id, '580')
		rej_chans = {'A12'};
	elseif strcmp(subj_id, '607') % check early and last large bursts
		rej_chans = {'A28', 'D10'};
	elseif strcmp(subj_id, '621')
		rej_chans = {'B32'};
	elseif strcmp(subj_id, '631')
		rej_chans = {'D28'};
	else
		rej_chans = {};
	end
	rej_chan_inds = [];
	for rej_chan = rej_chans
		labs = {EEG.chanlocs.labels};
		rej_chan_inds(end+1) = find(strcmp(labs,rej_chan));
	end
	EEG = pop_interp(EEG, rej_chan_inds, 'spherical');
	EEG.etc.pipeline{end+1} =  ['Channels  removed and interpolated: ', num2str(rej_chan_inds)];

	% copy ic results
	EEG.icaweights = ic_EEG.icaweights;
	EEG.icasphere  = ic_EEG.icasphere;
	EEG = eeg_checkset(EEG, 'ica');
	EEG.etc.pipeline{end+1} =  ['ICs copied from: ', ic_EEG.setname];


	% linked mastoid reference
	EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','off');
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'Linked-mastoid reref, M1 and M2 removed';

	% align head model, warp to fiducials and Cz(45)->C21(85), Iz(88)->D23(119)
	EEG = headfit(EEG,subj_id);

	EEG.setname = [old_setname,'_rej'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end