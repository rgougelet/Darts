clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'/eeglab2019_0/'])
eeglab;
close all;
data_dir = [script_dir,'/data/'];
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};
new_srate = 32;

% preprocess sets
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)

	% load dataset
  subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg.set'];
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% exclude channels on arm
	EEG = pop_select( EEG,'nochannel', {'EXT7', 'EXT8', 'EXG7', 'EXG8'});
	
	% optimize head center
	EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
	
	% hp filter
	EEG = pop_eegfiltnew(EEG, 1, 0, 1650, 0, [], 0);

	% notch filter
	EEG = pop_eegfiltnew(EEG, 59.5,60.5, [],true);
	
	% resample
	if EEG.srate ~= new_srate % keep old srate if equivalent
		EEG = pop_resample( EEG, new_srate, 0.8, 0.4);
	end
	
	% save set
	EEG.setname = [old_setname,'_',num2str(new_srate)];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
