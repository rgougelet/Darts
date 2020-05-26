clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'eeglab/plugins/firfilt2.4/'])
addpath(genpath([script_dir,'deps/']))
data_dir = [script_dir,'data/'];
addpath(genpath(data_dir))

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
baseline_length_in_samples = 128;
offset_in_samples = 384; % cut off end of trial to account for throw artifacts

each_clusters_label = {'Front','Back'};
labels_of_channels_in_each_cluster = {...
{'C11','C12','C6','C10', 'C15', 'C16','C17'},...
{'D7','D3','D6','D11','D12','D13', 'D8'}};