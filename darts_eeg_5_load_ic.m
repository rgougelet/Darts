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