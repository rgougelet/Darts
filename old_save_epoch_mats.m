rmpath('/data/common/matlab/eeglab')
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
addpath([script_dir,'eeglab/']);
addpath(genpath([script_dir,'deps/']));
% eeglab
%%
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
	} ;
data_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/data/';
pipe_dir = 'IIR HP 1 Hz Pass Edge - Notch Filters/';

for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir, pipe_dir,subj_id,'*_512.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',subj_set.folder);
	
	% load trial timings
	load([data_dir,'Latencies/', subj_id,'_eeg_',num2str(EEG.srate),'_latencies']);

	% get labels and indices for channels
	front_chan_i = [];
	for chan_lab = {'C11','C12','C6','C10', 'C15', 'C16','C17'}
		front_chan_i = [front_chan_i, find(strcmpi(chan_lab,{EEG.chanlocs.labels}))];
	end
	back_chan_i = [];
	for chan_lab = {'D7','D3','D6','D11','D12','D13', 'D8'}
		back_chan_i = [back_chan_i, find(strcmpi(chan_lab,{EEG.chanlocs.labels}))];
	end
	uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
	lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
	lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
	rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
	
	% get cluster means and veog
	prefront = mean(EEG.data(front_chan_i,:),1);
	preback = mean(EEG.data(back_chan_i,:),1);
	uveog = EEG.data(uveog_i,:);
	lveog = EEG.data(lveog_i,:);
	lheog = EEG.data(lheog_i,:);
	rheog = EEG.data(rheog_i,:);
	veog = uveog-lveog;
	heog = lheog-rheog;
	
	epochs = {};
for event_i = 1:length(start_event_latencies)
		epoch = [];
		epoch.inds = start_event_latencies(event_i):end_event_latencies(event_i)-384;% 384 to correct for motion artifacts
		epoch.prefront = prefront(:,epoch.inds)-mean(prefront(:,epoch.inds),2); 
		epoch.preback = preback(:,epoch.inds)-mean(preback(:,epoch.inds),2); 
		epoch.veog = veog(:,epoch.inds)-mean(veog(:,epoch.inds),2); 
		epoch.heog = heog(:,epoch.inds)-mean(heog(:,epoch.inds),2); 

		% eog regression
		Y = epoch.prefront';
		X = [ones(size(epoch.prefront)); epoch.veog; epoch.heog]';
		X_inv = pinv(X);
		b_front = X_inv*Y;
		epoch.postfront = (epoch.prefront'-X*b_front)';
		

		Y = epoch.preback';
		X = [ones(size(epoch.preback)); epoch.veog; epoch.heog]';
		X_inv = pinv(X);
		b_back = X_inv*Y;
		epoch.postback = (epoch.preback'-X*b_back)';


		epochs{event_i} = [epoch.postfront;epoch.postback];
	end
	save([data_dir,subj_id],'epochs')
end