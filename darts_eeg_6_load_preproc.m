clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath(genpath([script_dir,'deps/']))
data_dir = [script_dir,'data/Last Full Pipeline/'];
addpath(genpath(data_dir))

% eeglab nogui;
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

cases = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 0 0; 1 0 1; 1 1 1; 1 1 0];
for case_i = 1:length(cases(1:4,:))
	is_fir = cases(case_i,1);
	reject_labeled_ics = cases(case_i,2);
	use_laplacian = cases(case_i,3);
	if is_fir
		out_file = ['FIR -'];
		if reject_labeled_ics
			out_file = [out_file,' Reject ICs -'];
			if use_laplacian
				out_file = [out_file,' Laplacian'];
			else
				out_file = [out_file,' Linked Mastoid'];
			end
		else
			out_file = [out_file,' Keep ICs -'];
			if use_laplacian
				out_file = [out_file,' Laplacian'];
			else
				out_file = [out_file,' Linked Mastoid'];
			end
		end
	else
		out_file = ['IIR -'];
		if reject_labeled_ics
			out_file = [out_file,' Reject ICs -'];
			if use_laplacian
				out_file = [out_file,' Laplacian'];
			else
				out_file = [out_file,' Linked Mastoid'];
			end
		else
			out_file = [out_file,' Keep ICs -'];
			if use_laplacian
				out_file = [out_file,' Laplacian'];
			else
				out_file = [out_file,' Linked Mastoid'];
			end
		end
	end
	out_dir = [script_dir,'data/',out_file,'/'];
	mkdir(out_dir);
	out_file = [out_dir,out_file, '.mat'];
	addpath(genpath(out_dir))
	
	
	%% apply preprocessing to pre-ica data
	% parfor compatible
	parfor subj_i = 1:length(subjs_to_include)
		eeglab nogui;
		
		% load post-ica dataset
		subj_id = subjs_to_include{subj_i};
		subj_set = dir([data_dir, subj_id,'*_ic.set']);
		icEEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
		
		% load pre-ica dataset
		subj_set = dir([data_dir, subj_id,'*_',num2str(srate),'.set']);
		EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
		old_setname = EEG.setname;
		
		% remove rejected channels and interp
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
		EEG.etc.pipeline{end+1} =  ['Channel inds removed and interpolated: ', num2str(rej_chan_inds)];
		
		% avg ref
		EEG.nbchan = EEG.nbchan+1;
		EEG.data(end+1,:) = zeros(1, EEG.trials*EEG.pnts);
		EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
		EEG = pop_reref(EEG, []);
		EEG = pop_select( EEG,'nochannel',{'initialReference'});
		EEG.etc.pipeline{end+1} =  'Average reref.';
		
		% copy ic results
		EEG.icaweights = icEEG.icaweights;
		EEG.icasphere  = icEEG.icasphere;
		EEG = eeg_checkset(EEG, 'ica');
		EEG.etc.pipeline{end+1} =  ['ICs copied from: ', icEEG.setname];
		
		% redo eog pipeline
		uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
		lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
		lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
		rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
		
		in = load([data_dir, subj_id,'_eog']);
		mx = in.mx;	sdx = in.sdx;	B_inv = in.B_inv;
		Y = [EEG.data(uveog_i,:); EEG.data(lveog_i,:);EEG.data(lheog_i,:);EEG.data(rheog_i,:)]; % works best
		X = EEG.icaact(:,:);
		X = (X-mx)./sdx; % normalize
		X_hat = B_inv*Y; % back project from eog
		new_icaact= sdx.*(X-X_hat)+mx; %denormalize
		new_data = EEG.icawinv*new_icaact; % forward project to channels
		EEG.icaact = reshape(new_icaact,size(EEG.icaact));
		EEG.data = reshape(new_data, size(EEG.data));
		EEG.etc.pipeline{end+1} = 'ICs reweighted and EOG rejected';
		EEG.etc.eog.B_inv = B_inv;
		EEG.etc.eog.sdx = sdx;
		EEG.etc.eog.mx = mx;
		
		% linked-mastoid rereference
		EEG = pop_reref(EEG, {'M1','M2'}, 'keepref','off');
		EEG = eeg_checkset(EEG);
		EEG.etc.pipeline{end+1} =  'Linked-mastoid reref, M1 and M2 removed';
		EEG = pop_select( EEG,'nochannel', {'UVEOG', 'LVEOG', 'LHEOG', 'RHEOG'});
		EEG = eeg_checkset(EEG);
		EEG.etc.pipeline{end+1} =  'EOG channels removed';
		
		% copy over dipfit and label data
		subj_set = dir([data_dir, subj_id,'*_lab.set']);
		labEEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
		EEG.dipfit = labEEG.dipfit;
		EEG.etc.ic_classification.ICLabel = labEEG.etc.ic_classification.ICLabel;
		EEG = eeg_checkset(EEG);
		
		if reject_labeled_ics
			lab = labEEG.etc.ic_classification.ICLabel; %#ok<*UNRCH>
			classif = [(1:length(lab.classifications))', lab.classifications];
			rej_i = classif(:,2)<.75;
			% 	classif(rej_i,:) = [];
			% 	figure; heatmap(sortrows(classif(:,2:end),2))
			EEG = pop_subcomp( EEG, find(rej_i), 0);
			EEG.etc.pipeline{end+1} =  'Non-brain ICs removed';
			EEG.etc.pipeline{end+1} =  find(rej_i);
		end
		if use_laplacian
			%	apply laplacian
			X = [EEG.chanlocs.X];
			Y = [EEG.chanlocs.Y];
			Z = [EEG.chanlocs.Z];
			EEG.data = laplacian_perrinX(EEG.data,X,Y,Z);
			EEG.etc.pipeline{end+1} =  'Re-referenced to Laplacian';
			% save channel data
			EEG.setname = [old_setname,'_lap_6'];
			EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',out_dir,' at ', datestr(now)];
			pop_saveset(EEG, 'filename', EEG.setname,'filepath', out_dir);
		else
			% save channel data
			EEG.setname = [old_setname,'_lm_6'];
			EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',out_dir,' at ', datestr(now)];
			pop_saveset(EEG, 'filename', EEG.setname,'filepath', out_dir);
		end
		
		%% reduce to clusters
		chan_clusts = {...
			{'C11','C12','C6','C10', 'C15', 'C16','C17'},...% front middle
			{'D7','D3','D6','D11','D12','D13', 'D8'}... % back middle
			};
		clust_labs = {'Front Middle',  'Back Middle'};
		
		oldEEG = EEG;
		newEEG = eeg_emptyset();
		newEEG.times  = oldEEG.times;
		newEEG.srate  = srate;
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
			
			% add clusters as channels
			newEEG.data(end+1,:)  = mean(oldEEG.data(chan_i,:),1);
		end
		newEEG = eeg_checkset(newEEG);
		newEEG.etc.pipeline{end+1} =  'Reduced to clusters';
		newEEG.etc.pipeline{end+1} =  chan_clusts;
		
		% save set
		if use_laplacian
			EEG = newEEG;
			EEG.setname = [old_setname,'_lap_6_clust'];
			EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',out_dir,' at ', datestr(now)];
			pop_saveset(EEG, 'filename', EEG.setname,'filepath', out_dir);
		else
			EEG = newEEG;
			EEG.setname = [old_setname,'_lm_6_clust'];
			EEG.etc.pipeline{end+1} =  ['Saved as ',EEG.setname,' to ',out_dir,' at ', datestr(now)];
			pop_saveset(EEG, 'filename', EEG.setname,'filepath', out_dir);
		end
	end
end