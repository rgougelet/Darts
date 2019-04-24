clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis\';
cd(script_dir);
addpath('.\eeglab13_6_5b')
data_dir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis\data\';
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};

eeglab;
close all;
target_absent = [];
target_present = [];
srate = 64;
epoch_half_width = 3; % value of 3 would be -3 pre-event to 3 post-event, must be symmetric

for subj_i = 1:length(subjs_to_include)
	start = tic;
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir, subj_id,'*_laplace_epochs.set']);
	
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	target_absent_types = {};
	target_present_types = {};
	for event_i = 1:length(EEG.event)
		event_type = EEG.event(event_i).type;
		if length(event_type) <= 3
			continue
		end
		if strcmp(event_type(4),'0')
			target_absent_types{end+1} = event_type;
		end
		if strcmp(event_type(4),'1')
			target_present_types{end+1} = event_type;
		end
	end
	
	EEG_target_absent = pop_epoch( EEG, target_absent_types, [-epoch_half_width  epoch_half_width], 'newname', EEG.setname, 'epochinfo', 'yes');
	target_absent = cat(3,target_absent, EEG_target_absent.data);
	
	EEG_target_present = pop_epoch( EEG, target_present_types, [-epoch_half_width  epoch_half_width], 'newname', EEG.setname, 'epochinfo', 'yes');
	target_present = cat(3,target_present, EEG_target_present.data);

end

parsave('trial_data',{target_absent,target_present}, {'target_absent','target_present'})

% trash
	
% 	chan_i = 3;
% 	% Audio Active Rare-Frequent
% 	[aa_ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
% 		newtimef( {EEG_audio_active_rare.data(chan_i,:,:),EEG_audio_active_frequent.data(chan_i,:,:)},...
% 		1500, [-1000  2000], EEG.srate, [1 0.5],...
% 		'freqs', 1:0.5:80,...
% 		'freqscale', 'log',...
% 		'baseline', [-500, 0],...
% 		'timesout', 400);
% 	close;
% 	% Audio Resting Rare-Frequent
% 	[ar_ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
% 		newtimef( {EEG_audio_resting_rare.data(chan_i,:,:),EEG_audio_resting_frequent.data(chan_i,:,:)},...
% 		1500, [-1000  2000], EEG.srate, [1 0.5],...
% 		'freqs', 1:0.5:80,...
% 		'freqscale', 'log',...
% 		'baseline', [-500, 0],...
% 		'timesout', 400);
% 	close;
% 	% Visual Active Rare-Frequent
% 	[va_ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
% 		newtimef( {EEG_visual_active_rare.data(chan_i,:,:),EEG_visual_active_frequent.data(chan_i,:,:)},...
% 		1500, [-1000  2000], EEG.srate, [1 0.5],...
% 		'freqs', 1:0.5:80,...
% 		'freqscale', 'log',...
% 		'baseline', [-500, 0],...
% 		'timesout', 400);
% 	close;
% 	% Visual Resting Rare-Frequent
% 	[vr_ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
% 		newtimef( {EEG_visual_resting_rare.data(chan_i,:,:),EEG_visual_resting_frequent.data(chan_i,:,:)},...
% 		1500, [-1000  2000], EEG.srate, [1 0.5],...
% 		'freqs', 1:0.5:80,...
% 		'freqscale', 'log',...
% 		'baseline', [-500, 0],...
% 		'timesout', 400);
% 	close;
% 	fig = figure('Position', [100 100 1024 786]);
% 	subplot(2,2,1); tftopo(aa_ersp{3},times,freqs); title('Pz; Audio Walking Rare>Frequent');
% 	subplot(2,2,3); tftopo(ar_ersp{3},times,freqs); title('Pz; Audio Sitting Rare>Frequent');
% 	subplot(2,2,2); tftopo(va_ersp{3},times,freqs); title('Pz; Visual Walking Rare>Frequent');
% 	subplot(2,2,4); tftopo(vr_ersp{3},times,freqs); title('Pz; Visual Sitting Rare>Frequent');
% 	saveas(fig, num2str(subj_i), 'tiffn')
% 	close;

