clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis\';
cd(script_dir);
addpath('.\eeglab13_6_5b')
data_dir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis\data\';
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};

eeglab;
close all;
for subj_i = 1:length(subjs_to_include)

	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_64_laplace.set'];
	
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	
	load([data_dir, subj_id,'_eeg_64_latencys'])
	
end