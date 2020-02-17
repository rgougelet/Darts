clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
% script_dir = 'G:/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
eeglab nogui;
data_dir = [script_dir,'data/'];
addpath(data_dir)

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
new_srate = 512;

% preprocess sets
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'_eeg.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;
	EEG.etc.pipeline =  {};
	
	% exclude channels on arm
	EEG = pop_select( EEG,'nochannel', {'EXT7', 'EXT8', 'EXG7', 'EXG8'});
	EEG.etc.pipeline{end+1} =  'Arm channels removed';
	
	% fix channel labels, results verified manually
	if strcmp(subj_id,{'627'})
		ex = 'EXG';
	else
		ex = 'EXT';
	end
	m1_i = find(strcmp({EEG.chanlocs.labels},[ex,'1']));
	m2_i = find(strcmp({EEG.chanlocs.labels},[ex,'2']));
	uveog_i= find(strcmp({EEG.chanlocs.labels},[ex,'3']));
	lveog_i = find(strcmp({EEG.chanlocs.labels},[ex,'4']));
	lheog_i = find(strcmp({EEG.chanlocs.labels},[ex,'5']));
	rheog_i = find(strcmp({EEG.chanlocs.labels},[ex,'6']))
	EEG=pop_chanedit(EEG, 'changefield',{m1_i,'labels','M1'});
	EEG=pop_chanedit(EEG, 'changefield',{m2_i,'labels','M2'});
	EEG=pop_chanedit(EEG, 'changefield',{uveog_i,'labels','UVEOG'});
	EEG=pop_chanedit(EEG, 'changefield',{lveog_i,'labels','LVEOG'});
	EEG=pop_chanedit(EEG, 'changefield',{lheog_i,'labels','LHEOG'});
	EEG=pop_chanedit(EEG, 'changefield',{rheog_i,'labels','RHEOG'});
	EEG.etc.pipeline{end+1} =  'Externals relabeled';
	EEG.etc.pipeline = EEG.etc.pipeline'; % make cell array easier to display

	% optimize head center
	EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
	EEG.etc.pipeline{end+1} =  'Head center optimized';
	
	% linked-mastoid reref
	EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','on');
	EEG.etc.pipeline{end+1} =  'Linked-mastoid reref';
	
	% highpass filter the data
	[b,a] = butter(5,1/(EEG.srate/2),'high');
	EEG.data = single(filtfilt(b,a,double(EEG.data')))';
	EEG.etc.pipeline{end+1} = ...
		['Butterworth HP, b:',num2str(b),' a:',num2str(a)];

	% notch filter the data
	for harm = 60:60:(EEG.srate/2)
		nyq = EEG.srate/2;
		bw = pi*(.5/nyq); % hz width
		s = 2; % half-amplitude (-6 dB) cutoff, 1/4 amp cutoff, s = 4
		zz = w/2; 
		d = sqrt(s^2-1)*zz;
		r = 1-d;
		t = harm*pi/nyq;
		b = [1 -2*cos(t) 1];
		a = [1 -2*r*cos(t) r^2];
		freqz(b,a, EEG.srate*20,EEG.srate)
		EEG.data = single(filtfilt(b,a,double(EEG.data(:,:)')))';
		EEG.etc.pipeline{end+1} = ...
			['IIR notch, b:',num2str(b),' a:',num2str(a)];
	end
	
	% resample
	if EEG.srate ~= new_srate % keep old srate if equivalent
		EEG = pop_resample( EEG, new_srate, 0.8, 0.4);
		EEG.etc.pipeline{end+1} =  ...
			['Resampled from ',num2str(EEG.srate),' to ',num2str(new_srate)];
	end
	
	% save set
	EEG.setname = [old_setname,'_',num2str(EEG.srate)];
	EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',data_dir,' at ', datestr(now)];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
