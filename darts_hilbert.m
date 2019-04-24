%% Script timing
start = datetime('now');
debugging = 1;

%% Initialize paths
warning('off','all')
scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
cd(scriptdir)
rmpath /data/common/matlab/eeglab/
addpath(genpath('./eeglab14_1_2b'))

if debugging == 1
	datapath = './filtered_32/';
	cd(datapath)
	sets = dir('*32.set');
	cd(scriptdir)
	outpath = './hilbert_ics_32/';
	mkdir(outpath)
else
	datapath = './filtered_512/';
	cd(datapath)
	sets = dir('*512.set');
	cd(scriptdir)
	outpath = './hilbert_ics_512/';
	mkdir(outpath)
end
warning('on','all')

%% For each subject:freq
parfor set_i = 1:length(sets)
	 
	% import component data
	eeg1 = pop_loadset('filename', sets(set_i).name, 'filepath', datapath);
	eeg1 = eeg_checkset( eeg1 );
	eeg1.icaact = (eeg1.icaweights*eeg1.icasphere)*eeg1.data(eeg1.icachansind,:);
	
	% import trial start and end indices
	try
		trial_time_mat =load([scriptdir,'/trial_inds_32/',eeg1.setname(1:20),'_AD_32.mat']);
	catch
		error('Could not find trial indices')
	end
	end_event_time_inds = trial_time_mat.end_event_time_inds;
	start_event_time_inds = trial_time_mat.start_event_time_inds;
	
	% compute entropy for components
	comp_entropys = nan(length(eeg1.icawinv),length(end_event_time_inds));
	for comp_i = 1:length(eeg1.icawinv)
		prctdone(comp_i-1,length(eeg1.icawinv))
		comp = eeg1.icaact(comp_i,:);
		hilsig = hilbert(comp);
		amp = abs(hilsig);
		
		%% Localize to only within trials. Signal length not an issue within components.
		for trial_ind = 1:length(end_event_time_inds)
			trial_data = double(amp(start_event_time_inds(trial_ind):(end_event_time_inds(trial_ind)-eeg1.srate/2)));
			trial_entropy = entropy(pspectrum(nmlz(trial_data),eeg1.srate));
			comp_entropys(comp_i,trial_ind) = trial_entropy;
		end
	end
	
	%% Identify most informative comps
	[min_entropys,min_entropy_inds] = mink(mean(comp_entropys,2),1);
	comp = eeg1.icaact(min_entropy_inds,:);
	hilsig = hilbert(comp);
	amp = abs(hilsig);
	phases = angle(hilsig);
	phase_c = cos(phases);
	phase_s = sin(phases);
	freq = eeg1.srate*diff(unwrap(phases))/(2*pi);
	spec_data = [amp;phase_c;phase_s;[0, freq]];
	comp_entropy=comp_entropys(min_entropy_inds,:);
	matname = [outpath,eeg1.setname,'_ic',num2str(min_entropy_inds, '%03d'),'.mat'];
	parsave(matname,{spec_data,comp_entropy},{'spec_data','comp_entropy'});
end

%% Script timing
stop = datetime('now');
disp(['Start ',datestr(start)])
disp(['End ',datestr(stop)])
disp(['Runtime ', char(duration(stop-start))])