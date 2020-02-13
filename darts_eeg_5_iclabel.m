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

for subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir, subj_id,'*_ch.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_EEG = EEG;

	% apply laplacian 
% 	X = [EEG.chanlocs.X];
% 	Y = [EEG.chanlocs.Y];
% 	Z = [EEG.chanlocs.Z];
% 	EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);

	% reduce to 4 channels
	% 	chans2keep = {'C11','C26','D7','D17'}; % Fz, Cz, Pz, Oz
	% 	EEG = pop_select(EEG,'channel',chans2keep);

	EEG = pop_iclabel(EEG, 'default');
	lab = EEG.etc.ic_classification.ICLabel;
	classes = lab.classes;
	[sort_c, sort_i] = sort(lab.classifications(:,1),'ascend');

	% save set
	EEG.setname = [EEG.setname,'_lab'];
	pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
