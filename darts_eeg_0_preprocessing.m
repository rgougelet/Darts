% how the data would be converted to lower sample rate
% not necessary if data shared already downsampled
clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/darts';
script_dir = 'C:\Users\Rob\Desktop\darts';
cd(script_dir);
rmpath('/data/common/matlab/eeglab')
addpath('./eeglab13_6_5b')
data_dir = './data/';
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};

% user input
new_srate = 32;

for subj_i = 1:length(subjs_to_include)
	start = tic;
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg.set'];
	
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	
	% load dataset
	EEG = eeg_checkset( EEG );
	old_setname = EEG.setname;
	
	% exclude channels on arm
	EEG = pop_select( EEG,'nochannel', {'EXT7', 'EXT8', 'EXG7', 'EXG8'});
	
% 	% high pass filter the data
% 	filt_freq = 0.6;
% 	cutoff_dist = 1;
% 	window_type = 'blackman';
% 	filt_ord = pop_firwsord(window_type, EEG.srate, cutoff_dist);
% 	EEG = pop_firws(EEG, 'fcutoff', filt_freq/EEG.srate, 'ftype', 'highpass', 'wtype', window_type, 'forder', filt_ord);
 	
	% resample
	if EEG.srate ~= new_srate % keep old srate if equivalent
		EEG = pop_resample( EEG, new_srate, 0.8, 0.4);
	end
	
	% apply cleanline
	EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan ,'computepower',0,'linefreqs', 60:60:(EEG.srate/2) ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',4);
	
	% optimize head center
	EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
	
	% save set
	EEG.setname = [setname_prefix,'_',num2str(new_srate)];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
