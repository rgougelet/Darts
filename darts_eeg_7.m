clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
% script_dir = 'G:/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab/')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'deps/'])
eeglab; close;
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

%%
for subj_i = 1:length(subjs_to_include)
	
	%% load post-ica dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir, subj_id,'*_ic.set']);
	icEEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	%% load pre-ica data
	subj_set = dir([data_dir, subj_id,'*_',num2str(srate),'.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_setname = EEG.setname;
	
	%% copy ic results
	EEG.icaweights = icEEG.icaweights;
	EEG.icasphere  = icEEG.icasphere;
	EEG = eeg_checkset(EEG, 'ica');
	EEG.etc.pipeline{end+1} =  ['ICs copied from: ', icEEG.setname];
	
	%% remove rejected channels and interp
	if strcmp(subj_id, '580')
		rej_chans = {'A12'};
	elseif strcmp(subj_id, '607')
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
	
	%% redo eog pipeline
	uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
	lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
	lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
	rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
	
	in = load([data_dir, subj_id,'_eog']);
	mx = in.mx;
	sdx = in.sdx;
	B_inv = in.B_inv;
	Y = [EEG.data(uveog_i,:); EEG.data(lveog_i,:);EEG.data(lheog_i,:);EEG.data(rheog_i,:)]; % works best
	X = EEG.icaact(:,:);
	X = (X-mx)./sdx; % normalize
	X_hat = B_inv*Y; % back project from eog
	new_icaact= sdx.*(X-X_hat)+mx; %denormalize
	new_data = EEG.icawinv*new_icaact; % forward project to channels
	EEG.icaact = reshape(new_icaact,size(EEG.icaact));
	EEG.data = reshape(new_data, size(EEG.data));
	EEG.etc.pipeline{end+1} = 'ICs reweighted';
	EEG.etc.eog.B_inv = B_inv;
	EEG.etc.eog.sdx = sdx;
	EEG.etc.eog.mx = mx;
	
	EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','off');
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'Linked-mastoid reref, M1 and M2 removed';
	EEG = pop_select( EEG,'nochannel', {'UVEOG', 'LVEOG', 'LHEOG', 'RHEOG'});
	EEG = eeg_checkset(EEG);
	EEG.etc.pipeline{end+1} =  'EOG channels removed';
	
	%% copy over dipfit and label data
	subj_set = dir([data_dir, subj_id,'*_lab.set']);
	labEEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	EEG.dipfit = labEEG.dipfit;
	EEG.etc.ic_classification.ICLabel = labEEG.etc.ic_classification.ICLabel;
	EEG = eeg_checkset(EEG);
	
	lab = EEG.etc.ic_classification.ICLabel;
	classif = [(1:length(lab.classifications))', lab.classifications];
	rej_i = classif(:,2)<.75;
	% 	classif(rej_i,:) = [];
	% 	figure; heatmap(sortrows(classif(:,2:end),2))
	EEG = pop_subcomp( EEG, find(rej_i), 0);
	EEG.etc.pipeline{end+1} =  'Non-brain ICs removed';
	EEG.etc.pipeline{end+1} =  find(rej_i);

	%% reduce to clusters
	chan_clusts = {...
		{'C11','C12','C6','C10', 'C15', 'C16','C17'},...% front middle
		{'D7','D3','D6','D11','D12','D13', 'D8'}... % back middle
		};
	clust_labs = {'Front Middle',  'Back Middle'};
	
	oldEEG = EEG;
	newEEG = eeg_emptyset();
	newEEG.times  = oldEEG.times;
	newEEG.srate  = srate; % Rounded actual sampling rate. Note that the unit of the time must be in second.
	newEEG.nbchan = length(chan_clusts);
	newEEG.pnts   = size(oldEEG.data,2);
	newEEG.etc = oldEEG.etc;

	for chan_clust_i = 1:length(chan_clusts)
		chan_clust = chan_clusts{chan_clust_i};
		clust_lab = clust_labs{chan_clust_i};
		
		chan_i = [];
		for chan_lab_i = 1:length(chan_clust)
			chan_lab = chan_clust{chan_lab_i};
			chan_i = [chan_i, find(strcmpi(chan_lab,{oldEEG.chanlocs.labels}))];
		end
		
		% apply laplacian
		% 	X = [EEG.chanlocs.X];
		% 	Y = [EEG.chanlocs.Y];
		% 	Z = [EEG.chanlocs.Z];
		% 	EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);
		
		% add clusters as channels
		newEEG.data(end+1,:)  = mean(oldEEG.data(chan_i,:),1);
	end
	newEEG = eeg_checkset(newEEG);
	newEEG.etc.pipeline{end+1} =  'Reduced to clusters';
	newEEG.etc.pipeline{end+1} =  chan_clusts;
	
	% save set
	newEEG.setname = [old_setname,'_clust'];
	EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',data_dir,' at ', datestr(now)];
	pop_saveset(EEG, 'filename', newEEG.setname,'filepath', data_dir);
end
