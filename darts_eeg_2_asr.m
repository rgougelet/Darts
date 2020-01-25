close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
% rmpath('/data/common/matlab/eeglab')
% addpath([script_dir,'eeglab2019_0/'])

% addpath('/data/common/matlab/eeglab')
% rmpath([script_dir,'eeglab2019_0/'])

% if ~exist('EEG','var')
% 	eeglab;
% end
close all;
data_dir = [script_dir,'data/'];
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};
srate = 512;

for subj_i = 1:length(subjs_to_include)

	% load dataset
  subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_',num2str(srate),'.set'];
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% load event latencies
	load([data_dir, subj_id,'_eeg_',num2str(srate),'_latencys'])
	old_EEG = EEG;
	
	% cut out epochs
	new_EEG = old_EEG;
	new_EEG.data = [];
	for event_i = 1:length(start_event_latencys)
		epoch = old_EEG.data(:,start_event_latencys(event_i):end_event_latencys(event_i)-384);
		new_EEG.data = [new_EEG.data,epoch-mean(epoch,2)];
	end
	EEG = new_EEG;
	EEG.pnts = length(new_EEG.data);
	EEG.event = [];
	EEG = eeg_checkset(EEG);
% 	% clean
% 	block_size_in_secs = 1;
% 	block_size_in_samps = floor(EEG.srate*block_size_in_secs);
% 
% 	quant_size_in_samp = block_size_in_samps;
% 	tic
% 	n_blocks = floor(length(EEG.data)/block_size_in_samps);
% 	data_mask = zeros(size(EEG.data));
% 	block_mask = [];
% 	block_ps = [];
% 	for chan_i = 1:EEG.nbchan
% 		chan = EEG.data(chan_i,:);
% 		chan_quants = quantile(chan,0:(1/quant_size_in_samp):1);
% 		for block_i = 1:n_blocks
% 			block_start_i = (block_i-1)*block_size_in_samps+1;
% 			block_end_i = (block_i)*block_size_in_samps;
% 			block = chan(block_start_i:block_end_i);
% 			[h,p] = kstest2(block,chan_quants,0.05/n_blocks);
% 			data_mask(chan_i,block_start_i:block_end_i) = h;
% 			block_mask(chan_i,block_i) = h;
% 			block_ps(chan_i,block_i) = p;
% 		end
% 	end
% 	toc
% 	figure;heatmap(block_mask,'GridVisible','off')
% caxis([0 1])
% 	figure;plot(sum(double(block_ps>(0.001/n_blocks*EEG.nbchan)),1))
% good_data = EEG;
% good_data.event = [];
% good_data.data(logical(data_mask))=NaN;
% bad_data = EEG;
% bad_data.event = [];
% bad_data.data(~logical(data_mask))=NaN;
% vis_artifacts(good_data,bad_data, 'windowlength',40);

% asr
tic
EEG = clean_asr(EEG); % uses a modified version of clean_asr in this repo only
toc
EEG.setname = [old_setname,'_asr'];
EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);

% 	try
	% avg ref
% 	EEG.nbchan = EEG.nbchan+1;
% 	EEG.data(end+1,:) = zeros(1, EEG.pnts);
% 	EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
% 	EEG = pop_reref(EEG, []);
% 	EEG = pop_select( EEG,'nochannel',{'initialReference'});
 
	% amica
% 	mkdir([data_dir,'amicaResults/' EEG.setname,'/']);
% 	runamica15(EEG, 'num_chans', EEG.nbchan,...
% 			'numprocs',24,...
% 			'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
% 	EEG.etc.amica  = loadmodout15([data_dir,'amicaResults/' EEG.setname,'/']);
% 	EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
% 	EEG.icaweights = EEG.etc.amica.W;
% 	EEG.icasphere  = EEG.etc.amica.S;
% 	EEG = eeg_checkset(EEG, 'ica');
% 	
% 	EEG.setname = [old_setname,'_ica'];
% 	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
% 	catch
% 	end
end