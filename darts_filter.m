%% Script timing
start = datetime('now');
disp(['Start ',datestr(start)])
debugging = 1;

%% Initialize paths
scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
cd(scriptdir)
rmpath /data/common/matlab/eeglab/
addpath(genpath('./eeglab14_1_2b'))

if debugging == 1
	max_freq = 15;
	datapath = './eeg_32/';
	cd(datapath)
	sets = dir('*.set');
	cd(scriptdir)
	outpath = './filtered_32/';
	mkdir(outpath)
	addpath(outpath)
else
	max_freq = 59;
	datapath = './pre_asr_plus_ics/';
	addpath(datapath)
	cd(datapath)
	sets = dir('*.set');
	cd(scriptdir)
	outpath = './filtered_512/';
	mkdir(outpath)
	addpath(outpath)
end

%% Filter each subject
for i = 1:length(sets)
		eeg1 = pop_loadset('filename', sets(i).name, 'filepath', datapath);
		eeg1 = eeg_checkset( eeg1 );
		cutoff_dist = 0.2; % Hz, -6 dB
		window_type = 'blackman';
		filt_ord = pop_firwsord(window_type, eeg1.srate, 2*cutoff_dist);
		for filt_freq = 1:max_freq
			eeg2.setname = [eeg1.setname,'_f',num2str(filt_freq, '%03d')];
			out_filename = [eeg2.setname '.set'];
			if exist(out_filename, 'file') == 2
				disp('File exists, continuing...');
				continue
			end
			eeg2 = pop_firws(eeg1, 'fcutoff', [filt_freq-cutoff_dist, filt_freq+cutoff_dist], 'ftype', 'bandpass', 'wtype', window_type, 'forder', filt_ord); 			eeg2 = pop_saveset(eeg2, 'filename', out_filename, 'filepath', outpath);
		end
end

disp("Done.")