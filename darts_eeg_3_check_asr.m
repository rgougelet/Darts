close all; clc; clear;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab/')
addpath([script_dir,'eeglab2019_0/'])

eeglab;
close all;
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
for subj_i = 1:length(subjs_to_include)

	% load asr dataset
  subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'*_asr.set']);
	clean_EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	subj_set = dir([data_dir,subj_id,'*512.set']);
	old_EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	 vis_artifacts(clean_EEG,old_EEG, 'WindowLength',40);

end