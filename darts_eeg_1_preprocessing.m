clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\darts\';
cd(script_dir);
addpath('.\eeglab13_6_5b')
data_dir = '.\data\';
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};

eeglab;
close all;
for subj_i = 1:length(subjs_to_include)
	start = tic;
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg.set'];
	
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	
	% user input
	new_srate = 64;
	sets = dir([data_dir,'/*_eeg.set']);
	
	% load dataset
	eeg = pop_loadset('filename', sets(set_i).name, 'filepath', data_dir);
	eeg = eeg_checkset( eeg );
	old_setname = eeg.setname;
	
	% exclude channels on arm
	eeg = pop_select( eeg,'nochannel',{'EXT7', 'EXT8', 'EXG7', 'EXG8'});
	
	% high pass filter the data
	filt_freq = 0.6;
	cutoff_dist = 1;
	window_type = 'blackman';
	filt_ord = pop_firwsord(window_type, eeg.srate, cutoff_dist);
	eeg = pop_firws(eeg, 'fcutoff', filt_freq/eeg.srate, 'ftype', 'highpass', 'wtype', window_type, 'forder', filt_ord);
	
	% resample
	if eeg.srate ~= new_srate % keep old srate if equivalent
		eeg = pop_resample( eeg, new_srate, 0.8, 0.4);
	end
	
	% apply cleanline
	eeg = pop_cleanline(eeg, 'bandwidth',2,'chanlist',1:eeg.nbchan ,'computepower',0,'linefreqs', 60:60:(eeg.srate/2) ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',4);
	
	% optimize head center
	eeg = pop_chanedit(eeg, 'eval','chans = pop_chancenter( chans, [],[]);');
	
	% save set
	EEG.setname = [setname_prefix,'_64'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
	
end
