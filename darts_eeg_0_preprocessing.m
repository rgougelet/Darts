clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath(genpath([script_dir,'deps/']))
data_dir = [script_dir,'data/'];
addpath(genpath(data_dir))
out_dir = [data_dir,datestr(now,'mmm-dd-yyyy_HH-MM-SS'),'/'];
mkdir(out_dir);
addpath(genpath(out_dir))
eeglab nogui;
diary([out_dir,'Log.txt']);
RAII.diary = onCleanup(@() diary('off'));

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
	};
new_srate = 512;

% preprocess sets
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'_eeg.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;
	EEG.etc.pipeline =  {};
	
	% exclude channels on arm
	EEG = pop_select( EEG,'nochannel', {'EXT7', 'EXT8', 'EXG7', 'EXG8'});
	EEG.etc.pipeline{end+1} =  'Arm channels removed';
	
	% fix channel labels, results verified manually
	if strcmp(subj_id,{'627'})
		ex = 'EXG';
	else
		ex = 'EXT';
	end
	m1_i = find(strcmp({EEG.chanlocs.labels},[ex,'1']));
	m2_i = find(strcmp({EEG.chanlocs.labels},[ex,'2']));
	uveog_i= find(strcmp({EEG.chanlocs.labels},[ex,'3']));
	lveog_i = find(strcmp({EEG.chanlocs.labels},[ex,'4']));
	lheog_i = find(strcmp({EEG.chanlocs.labels},[ex,'5']));
	rheog_i = find(strcmp({EEG.chanlocs.labels},[ex,'6']));
	EEG=pop_chanedit(EEG, 'changefield',{m1_i,'labels','M1'});
	EEG=pop_chanedit(EEG, 'changefield',{m2_i,'labels','M2'});
	EEG=pop_chanedit(EEG, 'changefield',{uveog_i,'labels','UVEOG'});
	EEG=pop_chanedit(EEG, 'changefield',{lveog_i,'labels','LVEOG'});
	EEG=pop_chanedit(EEG, 'changefield',{lheog_i,'labels','LHEOG'});
	EEG=pop_chanedit(EEG, 'changefield',{rheog_i,'labels','RHEOG'});
	EEG.etc.pipeline{end+1} =  'Externals relabeled';
	EEG.etc.pipeline = EEG.etc.pipeline'; % make cell array easier to display

	% optimize head center
	EEG = pop_chanedit(EEG, 'eval','chans = pop_chancenter( chans, [],[]);');
	EEG.etc.pipeline{end+1} =  'Head center optimized';
	
	% linked-mastoid reref
	EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','on');
	EEG.etc.pipeline{end+1} =  'Linked-mastoid reref';
	figure; pwelch(EEG.data(:,:)',5000,20,[],EEG.srate,'onesided');
	title(['Prefilter - ' subj_id]);
	saveas(gcf,[out_dir,subj_id,'_Prefilter.jpg']);
	
	% highpass filter the data
	% some subjects won't be fixed with filter, need
	% to trim the data to epochs because they are moving
	% a lot outside of them
	[EEG.data, sr, wp, ws, rp, rs, n] = iirsos.hp(EEG.data(:,:),EEG.srate,1.5,2,10e-3,0);
	figure; pwelch(EEG.data(:,:)',5000,20,[],EEG.srate,'onesided');
	title(['After highpass - ' subj_id]);
	saveas(gcf,[out_dir,subj_id,'_Highpass.jpg']);

	EEG.etc.pipeline{end+1} = ...
			['Butterworth SOS HP: ', ...
			num2str(sr),...
			', ',num2str(wp),...
			', ',num2str(ws),...
			', ', num2str(rp),...
			', ', num2str(rs), ...
			', ', num2str(n)];
		
	% notch filter the data
	% some subjects won't be fixed with filter, need
	% to trim the data to epochs because they are moving
	% a lot outside of them
	for harm = 60:60:(EEG.srate/2)
			[EEG.data, sr, wp, ws, rp, rs, n] = iirsos.bs(EEG.data,EEG.srate,[harm-.25,harm+.25],[harm-.5,harm+.5],10e-3,0);
			EEG.etc.pipeline{end+1} = ...
			['Butterworth SOS Notch: ', ...
			num2str(sr),...
			', ',num2str(wp),...
			', ',num2str(ws),...
			', ', num2str(rp),...
			', ', num2str(rs), ...
			', ', num2str(n)];
	end
	figure; pwelch(EEG.data(:,:)',5000,20,[],512,'onesided');
	title(['After notches - ' subj_id]);
	saveas(gcf,[out_dir,subj_id,'_Notches.jpg']);

	% resample
	if EEG.srate ~= new_srate % keep old srate if equivalent
		EEG = pop_resample( EEG, new_srate, 0.8, 0.4);
		EEG.etc.pipeline{end+1} =  ...
			['Resampled from ',num2str(EEG.srate),' to ',num2str(new_srate)];
	end
	
	% save set
	EEG.setname = [old_setname,'_',num2str(EEG.srate)];
	EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',out_dir,' at ', datestr(now)];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', out_dir);
end
