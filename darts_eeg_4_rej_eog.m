clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
data_dir = [script_dir,'data/'];
addpath(data_dir)
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
eeglab nogui;

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
%% parfor compatible
for subj_i = 1:length(subjs_to_include)
	close all;
%%
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

	%%
	X2 = [veog;heog];
	X4 = [EEG.data(uveog_i,:); EEG.data(lveog_i,:);EEG.data(lheog_i,:);EEG.data(rheog_i,:)]; % works best

	X = EEG.icaact(:,:);
	ms = mean(X,2);
	sds = std(X,[],2);
	X = (X-ms)./sds;
	Y = X4;
	Y = (Y-mean(Y,2))./std(Y,[],2);
	B = Y*pinv(X);
	B_inv = pinv(B);
	X_hat = B_inv*Y;
	new_icaact = sds.*(X-X_hat)+ms;
	new_data = EEG.icawinv*new_icaact;
	EEG.icaact = reshape(new_icaact,size(EEG.icaact));
	EEG.data = reshape(new_data, size(EEG.data));
% 	eegplot(new_icaact,'data2', EEG.icaact(:,:), 'dispchans',32,'winlength',10)
% 	eegplot(new_data,'data2', EEG.data(:,:), 'dispchans',32,'winlength',10)

	EEG.etc.pipeline{end+1} = 'ICs reweighted';
	EEG.etc.pipeline{end+1} = B_inv;

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
% 	eeglab redraw
% 	EEG = headfit(EEG,subj_id);
% 	EEG.etc.pipeline{end+1} =  'Channels co-registered using headfit.m';
% 	EEG = pop_multifit(EEG, [1:size(EEG.icaact,1)] ,'threshold',100,'rmout','on','plotopt',{'normlen','on'});
% 	EEG.etc.pipeline{end+1} =  'Dipfit run';

	EEG.setname = [old_setname,'_ch'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
	
end