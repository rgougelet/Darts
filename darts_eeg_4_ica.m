clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
data_dir = [script_dir,'data/'];
addpath(data_dir)
rmpath([script_dir,'/eeglab2019_0/'])
addpath('/data/common/matlab/eeglab')
eeglab;

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

% binica is parfor compatible
parfor subj_i = 1:length(subjs_to_include)

	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'*_asr.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% avg ref
	EEG.nbchan = EEG.nbchan+1;
	EEG.data(end+1,:) = zeros(1, EEG.pnts);
	EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
	EEG = pop_reref(EEG, []);
	EEG = pop_select( EEG,'nochannel',{'initialReference'});
	EEG.etc.pipeline{end+1} =  'Average reref';
 
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
	EEG = pop_runica(EEG, 'icatype','binica','extended',1,'interupt', 'off');
	EEG.setname = [old_setname,'_binica'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
	
end