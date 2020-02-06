close all; clc; clear;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab/')
addpath([script_dir,'eeglab/'])
data_dir = [script_dir,'data/'];
addpath(data_dir)
eeglab nogui;

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

	% load dataset
  subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_',num2str(srate),'.set'];
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% load event latencies
	load([data_dir, subj_id,'_eeg_',num2str(srate),'_latencies'])
	old_EEG = EEG;
	
	% cut out epochs
	new_EEG = old_EEG;
	new_EEG.data = [];
	for event_i = 1:length(start_event_latencies)
		epoch = old_EEG.data(:,start_event_latencies(event_i):end_event_latencies(event_i)-384);
		new_EEG.data = [new_EEG.data,epoch-mean(epoch,2)];
	end
	EEG = new_EEG;
	EEG.pnts = length(new_EEG.data);
	EEG.event = [];
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'Non-trial data removed';

	% separate into chunks for easier noise removal
	chnk_length_in_secs = .25;
	chnk_length_in_samps = EEG.srate*chnk_length_in_secs;
	n_chnks = floor(EEG.pnts/chnk_length_in_samps);
	for chnk_i = 1:n_chnks
		EEG.event(end+1).type = 'split';
		EEG.event(end).latency = 1+(chnk_i-1)*chnk_length_in_samps;
	end
	EEG = pop_editeventvals(EEG,'latency',1);
	EEG = pop_epoch(EEG,{'split'},[0  chnk_length_in_secs],'epochinfo', 'yes');
	rej = 1:length({EEG.event.type});
	rej = rej(strcmp({EEG.event.type},'split'));
	EEG = pop_editeventvals(EEG,'delete',rej);
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'Separated into 250 ms chunks for easier cleaning';

	% save
	EEG.setname = [old_setname,'_trim'];
	EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',data_dir,' at ', datestr(now)];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);

end