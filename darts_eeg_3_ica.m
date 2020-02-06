clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
% script_dir = 'G:/darts/';
cd(script_dir);
data_dir = [script_dir,'data/'];
addpath(data_dir)
rmpath([script_dir,'/eeglab/'])
addpath('/data/common/matlab/eeglab')
eeglab nogui;

subjs_to_include = {
% 	'571'
		'579'
% 	'580'
% 	'607'
% 	'608'
% 	'616'
% 	'619'
% 	'621'
% 	'627'
% 	'631'
	};
srate = 512;
%%
% binica is parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'*_trim.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% reject epochs
	rej_ep_inds = find(EEG.reject.rejmanual);
	EEG.etc.pipeline{end+1} =  ['Epochs removed: ', num2str(rej_ep_inds)];
	EEG = pop_rejepoch( EEG, rej_ep_inds ,0);
	
	% reject channels
	if strcmp(subj_id, '580')
		rej_chans = {'A12'};
	elseif strcmp(subj_id, '607') % also has large bursts early and late in recording, hoping they will be removed
		rej_chans = {'A28', 'D10'};
	elseif strcmp(subj_id, '621')
		rej_chans = {'B32'};
	elseif strcmp(subj_id, '631')
		rej_chans = {'D28'};
	else
		rej_chans = {};
	end
	
	rej_chan_inds = [];
	for rej_chan = rej_chans
		labs = {EEG.chanlocs.labels};
		rej_chan_inds(end+1) = find(strcmp(labs,rej_chan));
	end
	EEG = pop_interp(EEG, rej_chan_inds, 'spherical');
	EEG.etc.pipeline{end+1} =  ['Channels  removed and interpolated: ', num2str(rej_chan_inds)];
	
	% avg ref
	EEG.nbchan = EEG.nbchan+1;
	EEG.data(end+1,:) = zeros(1, EEG.trials*EEG.pnts);
	EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
	EEG = pop_reref(EEG, []);
	EEG = pop_select( EEG,'nochannel',{'initialReference'});
	EEG.etc.pipeline{end+1} =  'Average reref.';
	
	% amica
	% out_dir = [data_dir,'amicaResults/' EEG.setname,'/'];
	% mkdir(out_dir);
	% runamica15(EEG, 'num_chans', EEG.nbchan,...
	% 'numprocs',24,...
	% 'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1, 'outdir',out_dir);
	% EEG.etc.amica  = loadmodout15([data_dir,'amicaResults/' EEG.setname,'/']);
	% EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
	% EEG.icaweights = EEG.etc.amica.W;
	% EEG.icasphere  = EEG.etc.amica.S;
	% EEG = eeg_checkset(EEG, 'ica');
	% EEG.setname = [old_setname,'_amica'];
	% EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
	
	% binica
	cd(data_dir);
	EEG = pop_runica(EEG, 'icatype','binica','extended',1,'interupt', 'off', 'filenum', str2num(subj_id));
	cd(script_dir);
	EEG.setname = [old_setname,'_ic'];
	EEG.etc.pipeline{end+1} =  'BinICA run.';
	EEG.etc.pipeline{end+1} =  EEG.icaweights;
	EEG.etc.pipeline{end+1} =  EEG.icasphere;
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
	
end