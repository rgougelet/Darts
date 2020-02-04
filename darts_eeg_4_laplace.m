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
srate = 512;

% apply laplacian spatial filter
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_',num2str(srate),'.set'];
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	
	% apply laplacian and reduce to 4 channels
	chans2keep = {'C11','C26','D7','D17'}; % Fz, Cz, Pz, Oz
	X = [EEG.chanlocs.X];
	Y = [EEG.chanlocs.Y];
	Z = [EEG.chanlocs.Z];
	EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);
	EEG = pop_select(EEG,'channel',chans2keep);
	
	% save set
	EEG.setname = [EEG.setname,'_laplace'];
	pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end
