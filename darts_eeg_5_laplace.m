clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab/')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'deps/'])
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
srate = 512;

% apply laplacian spatial filter
% parfor compatible
for subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir, subj_id,'*_trim.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_EEG = EEG;
	% apply laplacian 
	X = [EEG.chanlocs.X];
	Y = [EEG.chanlocs.Y];
	Z = [EEG.chanlocs.Z];
	EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);
% 	eegplot(EEG.data,'srate',EEG.srate,'data2',old_EEG.data./std(old_EEG.data(:,:),0,2))
	% reduce to 4 channels
	% 	chans2keep = {'C11','C26','D7','D17'}; % Fz, Cz, Pz, Oz
	% 	EEG = pop_select(EEG,'channel',chans2keep);
	new = EEG;
	new.data = new.data(:,:)./std(new.data(:,:),0,2);
	old = old_EEG;
	old.data = old.data(:,:)./std(old.data(:,:),0,2);
	vis_artifacts_rjg(new,old,'equalize_channel_scaling',true);
	% save set
	EEG.setname = [EEG.setname,'_laplace'];
	pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
