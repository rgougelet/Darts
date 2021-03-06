clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'deps/'])
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

%% use iclabel plugin to identify non-brain ics
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir, subj_id,'*_ch.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_EEG = EEG;

	% apply laplacian (experimental)
% 	X = [EEG.chanlocs.X];
% 	Y = [EEG.chanlocs.Y];
% 	Z = [EEG.chanlocs.Z];
% 	EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);

	% reduce to 4 channels (experimental)
	% 	chans2keep = {'C11','C26','D7','D17'}; % Fz, Cz, Pz, Oz
	% 	EEG = pop_select(EEG,'channel',chans2keep);

	EEG = pop_iclabel(EEG, 'default');
	lab = EEG.etc.ic_classification.ICLabel;
	classes = lab.classes;
	[sort_c, sort_i] = sort(lab.classifications(:,1),'ascend');
% 	figure; eegplot(EEG.icaact(sort_i,:), 'dispchans',32, 'limits', [0 5], 'winlength',5)
	EEG.etc.pipeline{end+1} =  'Labels saved';

	% save set
	EEG.setname = [EEG.setname,'_lab'];
	EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',data_dir,' at ', datestr(now)];
	pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
