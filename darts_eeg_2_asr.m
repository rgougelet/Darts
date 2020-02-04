close all; clc; clear;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab/')
addpath([script_dir,'eeglab/'])

% rmpath([script_dir,'eeglab2019_0/'])
% addpath('/data/common/matlab/eeglab')

eeglab;
close all;
data_dir = [script_dir,'data/'];
addpath(data_dir)

subjs_to_include = {
		'571'
% 		'579'
% 		'580'
% 		'607'
% 		'608'
% 		'616'
% 		'619'
% 		'621'
% 		'627'
% 		'631'
	};
srate = 512;
for subj_i = 1:length(subjs_to_include)

	% load dataset
  subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_',num2str(srate),'.set'];
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% load event latencies
	load([data_dir, subj_id,'_eeg_',num2str(srate),'_latencys'])
	old_EEG = EEG;
	
	% cut out epochs
	new_EEG = old_EEG;
	new_EEG.data = [];
	for event_i = 1:length(start_event_latencys)
		epoch = old_EEG.data(:,start_event_latencys(event_i):end_event_latencys(event_i)-384);
		new_EEG.data = [new_EEG.data,epoch-mean(epoch,2)];
	end
	EEG = new_EEG;
	EEG.pnts = length(new_EEG.data);
	EEG.event = [];
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'Non-trial data removed';

% 	eegplot(EEG.data(1:5:end,:))
		vis_artifacts_rjg(EEG.data, 'WindowLength',40)

	% asr
% 	EEG = clean_asr(EEG,[],[],[],[],[],[],[],true); % uses gpu
% 	EEG.setname = [old_setname,'_asr'];
% 	EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',data_dir,' at ', datestr(now)];
% 	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);

end