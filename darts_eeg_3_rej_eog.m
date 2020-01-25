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

for subj_i = 1%:length(subjs_to_include)

	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'*_binica.set'];
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	old_setname = EEG.setname;
	
	% linked mastoid reference
	EEG.nbchan = EEG.nbchan+1;
	EEG.chanlocs(1,EEG.nbchan).labels = 'linked_mastoid';
	m_one_i = find(strcmp({EEG.chanlocs.labels},'M1'));
	m_two_i = find(strcmp({EEG.chanlocs.labels},'M2'));
	EEG.data(end+1,:) = (EEG.data(m_one_i,:)+EEG.data(m_two_i,:))/2;
	EEG = pop_reref(EEG, EEG.nbchan);
	EEG = pop_select( EEG,'nochannel',{'linked_mastoid'});
	EEG.etc.pipeline{end+1,:} =  'Linked mastoid reref';
 
	% amica
% 	out_dir = [data_dir,'amicaResults/' EEG.setname,'/'];
% 	mkdir(out_dir);
% 	runamica15(EEG, 'num_chans', EEG.nbchan,...
% 			'numprocs',24,...
% 			'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1, 'outdir',out_dir);
% 	EEG.etc.amica  = loadmodout15([data_dir,'amicaResults/' EEG.setname,'/']);
% 	EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
% 	EEG.icaweights = EEG.etc.amica.W;
% 	EEG.icasphere  = EEG.etc.amica.S;
% 	EEG = eeg_checkset(EEG, 'ica');

	EEG = pop_runica(EEG, 'icatype','binica','extended',1,'interupt', 'off');
	
	EEG.setname = [old_setname,'_binica'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
	
end