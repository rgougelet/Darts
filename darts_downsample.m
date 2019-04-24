[scriptdir,~,~] = fileparts(mfilename('fullpath'));
cd(scriptdir)
datapath = [scriptdir,'/pre_asr_plus_ics/'];
cd(datapath)
sets = dir('*.set');
outpath = [scriptdir,'/downsampled/'];
mkdir(outpath)
for i = 1:length(sets)
	EEG = pop_loadset('filename', sets(i).name, 'filepath', datapath);
	EEG = eeg_checkset( EEG );
	EEG = pop_resample( EEG, 32);
	EEG.setname = [sets(i).name(1:23),'_32'];
	pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', outpath);
end
