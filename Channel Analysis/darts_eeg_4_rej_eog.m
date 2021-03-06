clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'deps/'])
data_dir = [script_dir,'data/'];
addpath(data_dir)

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
%% use eog channels to remove eye blinks
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	eeglab nogui; % each proc needs eeglab dependencies

	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'*_ic.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;

	% get eog inds
	uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
	lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
	lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
	rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
	veog = EEG.data(uveog_i,:)-EEG.data(lveog_i,:);
	heog = EEG.data(lheog_i,:)-EEG.data(rheog_i,:);

	% EEG.icaact = EEG.icaweights*EEG.data
	% EEG.data = EEG.icawinv*EEG.icaact
	Y2 = [veog;heog]; % two bipolar channels
	Y4 = [EEG.data(uveog_i,:); ...
		EEG.data(lveog_i,:);...
		EEG.data(lheog_i,:);...
		EEG.data(rheog_i,:)]; % four unipolar channels, works best
	
	% normalize eeg data
	X = EEG.icaact(:,:);
	mx = mean(X,2);	sdx = std(X,[],2);
	X = (X-mx)./sdx; 
	
	% normalize eog data
	Y = Y4; % use four unipolar here
	my = mean(Y,2);	sdy = std(Y,[],2);
	Y = (Y-my)./sdy;
	
	% regression
	B = Y*pinv(X);
	B_inv = pinv(B);
	X_hat = B_inv*Y; % back project from eog
	
	new_icaact= sdx.*(X-X_hat)+mx; % denormalize
% 	eegplot(new_icaact, 'dispchans',32,'winlength',10, 'spacing', 30)
% 	eegplot(EEG.icaact(:,:), 'dispchans',32,'winlength',10, 'spacing', 30)

	new_data = EEG.icawinv*new_icaact; % forward project to channels
% 	eegplot(new_data, 'dispchans',32,'winlength',10, 'spacing', 30)
% 	eegplot(EEG.data(:,:), 'dispchans',32,'winlength',10, 'spacing', 30)
	
	EEG.icaact = reshape(new_icaact,size(EEG.icaact));
	EEG.data = reshape(new_data, size(EEG.data));

	EEG.etc.pipeline{end+1} = 'ICs reweighted and EOG rejected';
	EEG.etc.eog.B_inv = B_inv;
	EEG.etc.eog.sdx = sdx;
	EEG.etc.eog.mx = mx;
	parsave([data_dir, subj_id,'_eog'],{B_inv, sdx,mx},{'B_inv', 'sdx','mx'});

% 	% isolate oscillatory ICs (does not work, hard to isolate)
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
	EEG = headfit(EEG,subj_id);
	EEG.etc.pipeline{end+1} =  'Channels co-registered using headfit.m';
	EEG = pop_multifit(EEG, [1:size(EEG.icaact,1)] ,'threshold',100,'rmout','on','plotopt',{'normlen','on'});
	EEG.etc.pipeline{end+1} =  'Dipfit run';
	dipfit = EEG.dipfit;
	parsave([data_dir,subj_id,'_dipfit'],dipfit, 'dipfit');

	% save
	EEG.setname = [old_setname,'_ch'];
	EEG = pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);
end