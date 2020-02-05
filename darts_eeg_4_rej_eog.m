clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
data_dir = [script_dir,'data/'];
addpath(data_dir)
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'/eeglab2019_0/'])
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
	subj_set = dir([data_dir,subj_id,'*_binica.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	% get eog inds
	uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
	lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
	lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
	rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
	
	% linked mastoid reference
	EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','on');
	veog = EEG.data(uveog_i,:)-EEG.data(lveog_i,:);
	heog = EEG.data(lheog_i,:)-EEG.data(rheog_i,:);
	ic_EEG = EEG;

	% reject eog
	veog_rej_inds = abs(corr(veog',EEG.icaact'))>0.08;
	heog_rej_inds = abs(corr(heog',EEG.icaact'))>0.08;
	rej_inds = heog_rej_inds | veog_rej_inds;
	
	% apply
	subj_set = dir([data_dir,subj_id,'_eeg_512.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;
	
	EEG.icaweights = ic_EEG.icaweights  ;
	EEG.icasphere  = ic_EEG.icasphere;
	EEG = eeg_checkset(EEG, 'ica');
	
	EEG.etc.pipeline{end+1,:} =  'External channels removed';

	EEG.setname = [old_setname,'_rej'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end