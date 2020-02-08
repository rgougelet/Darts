clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
data_dir = [script_dir,'data/'];
addpath(data_dir)
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
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
%%
for subj_i = 1:length(subjs_to_include)
	close all;

	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'*_ic.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;
	old_EEG = EEG;

	% get eog inds
	uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
	lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
	lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
	rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
	veog = EEG.data(uveog_i,:)-EEG.data(lveog_i,:);
	heog = EEG.data(lheog_i,:)-EEG.data(rheog_i,:);

	% EEG.icaweights*EEG.data = EEG.icaact
	% EEG.data = EEG.icawinv*EEG.icaact
	
	X2 = [veog;heog];
	X4 = [EEG.data(uveog_i,:); EEG.data(lveog_i,:);EEG.data(lheog_i,:);EEG.data(rheog_i,:)];

% 	% uses bipolar eog channels and channels
% 	Y = EEG.data(1:end-4,:);
% 	X = X2;
% 	B = Y*pinv(X);
% 	new_Y = Y-B*X;
% 	eegplot([Y;X4],'data2',[new_Y;X4] , 'dispchans', 32, 'winlength',10);

% 	% uses unipolar eog channels and channels
% 	Y = EEG.data(1:end-4,:);
% 	X = X4;
% 	B = Y*pinv(X);
% 	new_Y = Y-B*X;
% 	eegplot([Y;X4],'data2',[new_Y;X4] , 'dispchans', 32, 'winlength',10);

% 	% uses bipolar eog channels and ICs
% 	Y = EEG.icaact(:,:);
% 	X = X2;
% 	B = Y*pinv(X);
% 	new_Y = Y-B*X;
% 	eegplot(EEG.data(:,:),'data2',EEG.icawinv*(new_Y), 'dispchans', 32, 'winlength',10);

	% uses all unipolar eog channels and ICs, seems best
	Y = EEG.icaact(:,:);
	X = X4;
	B = Y*pinv(X);
	new_Y = Y-B*X;
	EEG.icaact = EEG.icaact-reshape(B*X,size(EEG.icaact));
% 	eegplot(EEG.data(:,:),'data2',EEG.icawinv*(new_Y), 'dispchans', 32, 'winlength',10);

	EEG.etc.pipeline{end+1} = 'ICs reweighted';
	EEG.etc.pipeline{end+1} = B;

% 	% isolate oscillatory ICs (failed)
% 	[icx,f] = pwelch(diff(EEG.icaact(:,:)',1),EEG.srate*100,EEG.srate*50,[],EEG.srate,'onesided');
% 	eegplot(diff(EEG.icaact(:,:),2));
% 	plot(f,icx); xlim([0 50])
% 	theta_i = f>=3 & f< 8;
% 	alpha_i = f>=8 & f< 12;
% 	beta_i = f>=12 & f< 30;
% 	logamma_i = f>=30 & f< 80;
% 	higamma_i = f>=80 & f< 256;
% 	alpha_ratios = mean(icx(theta_i,:),1)./mean(icx(~theta_i,:),1);%thetza
% 	[~, sortalph ] = sort(alpha_ratios,'descend');
% 	figure; plot(f,icx(:,sortalph(1:5))); xlim([0 50])	
% 	eegplot(EEG.icaact(sortalph,:));
 
  % linked mastoid rereference
	EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','off');
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'Linked-mastoid reref, M1 and M2 removed';

	% remove eog channels
	EEG = pop_select( EEG,'nochannel', {'UVEOG', 'LVEOG', 'LHEOG', 'RHEOG'});
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'EOG channels removed';

	% align head model, warp to fiducials and Cz(45)->C21(85), Iz(88)->D23(119)
	eeglab redraw
	EEG = headfit(EEG,subj_id);
	in_cell = {'mesh',...
		'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_vol.mat',...
		'transform',...
		[],...
		'chaninfo1',...
		EEG.chaninfo,...
		'helpmsg',...
		'on',...
		'transform',...
		EEG.dipfit.coord_transform;
		};
	ref_locs = 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\elec\standard_1020.elc';

	[~, EEG.dipfit.coord_transform] = coregister(EEG.chanlocs,ref_locs,in_cell{:});
	EEG.etc.pipeline{end+1} =  'Channels co-registered using headfit.m then manually';
	EEG.etc.pipeline{end+1} = EEG.dipfit.coord_transform;

	EEG.setname = [old_setname,'_ch'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
	continue
	%% load pre-ica data
	ic_EEG = EEG;
	subj_set = dir([data_dir,subj_id,'_eeg_512.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;

	% reject channels
	if strcmp(subj_id, '580')
		rej_chans = {'A12'};
	elseif strcmp(subj_id, '607') % check early and last large bursts
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

	% copy ic results
	EEG.icaweights = ic_EEG.icaweights;
	EEG.icasphere  = ic_EEG.icasphere;
	EEG = eeg_checkset(EEG, 'ica');
	EEG.etc.pipeline{end+1} =  ['ICs copied from: ', ic_EEG.setname];

	% rej eog ics
	EEG = pop_subcomp( EEG, eog_rej_inds, 0);
	EEG.etc.pipeline{end+1} =  ['ICs rejected: ', num2str(eog_rej_inds)];
	EEG = eeg_checkset(EEG, 'ica');

	% linked mastoid reference
	EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','off');
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'Linked-mastoid reref, M1 and M2 removed';

	% remove eog channels
	EEG = pop_select( EEG,'nochannel', {'UVEOG', 'LVEOG', 'LHEOG', 'RHEOG'});
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'EOG channels removed';

	% align head model, warp to fiducials and Cz(45)->C21(85), Iz(88)->D23(119)
	EEG = headfit(EEG,subj_id);
	EEG.etc.pipeline{end+1} =  'Channels co-registered using headfit.m';

	EEG.setname = [old_setname,'_rej'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end