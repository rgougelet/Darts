clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\darts\';
cd(script_dir);
addpath('.\eeglab13_6_5b')
data_dir = '.\data\';
addpath(data_dir)

% epoch_half_width = 3; % value of 3 would be -3 pre-event to 3 post-event, must be symmetric
subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};

eeglab;
close all;
for subj_i = 1:length(subjs_to_include)
	start = tic;
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_64.set'];
	
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	
% 	target_onset_types = {};
% 	for event_i = 1:length(EEG.event)
% 		event_type = EEG.event(event_i).type;
% 		if length(event_type) <= 3
% 			continue
% 		end
% 		if strcmp(event_type(end-2:end),'003')
% 			target_onset_types{end+1} = event_type;
% 		end
% 	end
	
	% near Fz
	chan1 = 'C11';
	chan1_neighbors = {'C6', 'C10', 'C15', 'C16', 'C17', 'C12'};
	
	% near Cz
	chan2 = 'C26';
	chan2_neighbors = {'C21', 'C25', 'C30', 'D1', 'C31', 'C27'};
	
	% near Pz
	chan3 = 'D7';
	chan3_neighbors = {'D3', 'D6', 'D11', 'D12', 'D13', 'D8'};
	
	% near Oz
	chan4 = 'D17';
	chan4_neighbors = {'D12', 'D16', 'D22', 'D23', 'D24', 'D18'};
	
	EEG_1_refs = pop_select( EEG,'channel',chan1_neighbors);
	EEG_2_refs = pop_select( EEG,'channel',chan2_neighbors);
	EEG_3_refs = pop_select( EEG,'channel',chan3_neighbors);
	EEG_4_refs = pop_select( EEG,'channel',chan4_neighbors);
	
	EEG_out = pop_select( EEG,'channel',{chan1, chan2, chan3, chan4});
	chan1idx = strcmpi(chan1,{EEG_out.chanlocs.labels});
	chan2idx = strcmpi(chan2,{EEG_out.chanlocs.labels});
	chan3idx = strcmpi(chan3,{EEG_out.chanlocs.labels});
	chan4idx = strcmpi(chan4,{EEG_out.chanlocs.labels});
	
	EEG_out.data(chan1idx,:) = EEG_out.data(chan1idx,:)-mean(EEG_1_refs.data)/length(chan1_neighbors);
	EEG_out.data(chan2idx,:) = EEG_out.data(chan2idx,:)-mean(EEG_2_refs.data)/length(chan2_neighbors);
	EEG_out.data(chan3idx,:) = EEG_out.data(chan3idx,:)-mean(EEG_3_refs.data)/length(chan3_neighbors);
	EEG_out.data(chan4idx,:) = EEG_out.data(chan4idx,:)-mean(EEG_4_refs.data)/length(chan4_neighbors);
	
	EEG = EEG_out;
	EEG.setname = [EEG.setname,'_laplace'];
% 	EEG = pop_epoch( EEG, target_onset_types, [-epoch_half_width  epoch_half_width], 'newname', EEG.setname, 'epochinfo', 'yes');
% 	EEG = pop_jointprob(EEG,1,1:4 ,5,5,0,0,0,[],0);
% 	EEG = pop_rejkurt(EEG,1,1:4 ,5,5,0,0,0,[],0);
% 	EEG = pop_rejepoch( EEG,find(EEG.reject.rejjp==1 | EEG.reject.rejkurt==1) ,0);
	
	pop_saveset(EEG, 'filename', EEG.setname,'filepath', data_dir);

end




