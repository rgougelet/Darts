clear; close all; clc;
script_dir = '/home/rgougelet/Desktop/darts/';
cd(script_dir);
rmpath('/data/common/matlab/eeglab')
addpath('./eeglab14_1_2b')
eeglab;
close all;
data_dir = './data/';
addpath(data_dir)

subjs_to_include = {'571','579', ...
	'607', '608', '616', '619', '621', '627', '631', '580'};

% user input
new_srate = 512;

for subj_i = 1:length(subjs_to_include)

	% load dataset
    subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg.set'];
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% exclude channels on arm
	EEG = pop_select( EEG,'nochannel', {'EXT7', 'EXT8', 'EXG7', 'EXG8'});
	
	% optimize head center
	EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
	
	% resample
	if EEG.srate ~= new_srate % keep old srate if equivalent
		EEG = pop_resample( EEG, new_srate, 0.8, 0.4);
	end
	
	% apply cleanline
% 	EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan ,'computepower',0,'linefreqs', 60:60:(EEG.srate/2) ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',4);
	
	% save set
	EEG.setname = [old_setname,'_',num2str(new_srate)];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
