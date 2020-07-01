clear
script_dir= '/data/common/mobi/Experiments/Darts/Analysis/darts/';
addpath(script_dir)
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'eeglab/plugins/firfilt2.4/'])
addpath(genpath([script_dir,'deps/']))
data_dir = [script_dir,'data/'];
addpath(genpath(data_dir))

srate = 512;
% labels_of_each_cluster = {'Front','Back'};
% labels_of_channels_in_each_cluster = {...
% {'C11','C12','C6','C10', 'C15', 'C16','C17'},...
% {'D7','D3','D6','D11','D12','D13', 'D8'}};

baseline_length_in_samples = 128;
offset_in_samples = 384; % cut off end of trial to account for throw artifacts

% pipe_name= 'IIR HP 1 Hz Pass Edge - Notch Filters/';
% pipe_name ='IIR BP 1 Hz Pass Edge - Lower Order - Notch Filters - Monopolar';
% pipe_name ='IIR BP 1 Hz Pass Edge - Lower Order - Notch Filters - Bipolar';
% pipe_name = 'IIR BP 1 Hz Pass Edge - Lower Order - Notch Filters - Laplacian';
dbstop if error
pipe_dir = [data_dir,pipe_name,'/'];

%% create subject structure

% size(pipeline.Subjects(1).Clusters(1).Data.EOG.Raw,1)
pipeline = create_pipeline_structure(data_dir,pipe_name,srate,1);
% get clusters % considering single vs bipolar channel eog
% pipeline = get_clusters(pipeline,labels_of_each_cluster,labels_of_channels_in_each_cluster,'monopolar');
% pipeline = get_clusters(pipeline,labels_of_each_cluster,labels_of_channels_in_each_cluster,'bipolar');
% get epochs
pipeline = get_epochs(pipeline,baseline_length_in_samples,offset_in_samples);
% get trial amps
pipeline = get_trial_amps(pipeline);
% get behavioral data
pipeline = get_behav_data(pipeline);
% save([pipe_dir,'pipeline_structure.mat'],'pipeline', '-v7.3')
% load([pipe_dir,'pipeline_structure.mat'])

% plot epochs
% plot_epochs(pipeline,1);
%%
cluster_stats(pipeline,2)