%% cleanEEG
% EEG data preprocessing should be fast and easy, require no special algorithms,
% and be able to reasonably recover data from typical EEG datasets.
% It strictly uses well-known summary statistics and fast visualizations
% to make cleaning the data intuitive, easy, and most importantly, quick.

clc
disp('Welcome to the EEG cleaning utility...')
main_str = '';
temp_EEG = [];

if exist('EEG', 'var')
	disp('Using EEG object already in memory')
end
if ~exist('EEG', 'var')
	error(['You must load a dataset into a structure named EEG. Try running',...
		' EEGLAB and loading a dataset'])
end
while isempty(main_str)
	if ~exist('gui','var') || ~isvalid(gui)
		gui = figure('Name','EEG Cleaning Utility','NumberTitle','off');
		set(gui,'WindowStyle','docked')
	else
		figure(gui)
	end
	if isempty(temp_EEG)
		temp_EEG = EEG;
	end
	prompt = [...
		'\nMain Menu\n',...
		' 1. Load new dataset (make sure to save first, if needed)\n'...
		' 2. Remove epochs/chunks by standard deviation\n'...
		' 3. Remove epochs/chunks by amplitude\n',...
		' 4. Quick plot EEG data\n',...
		' 5. Remove channels via pairwise correlation\n'...
		' 6. Plot data with vis_artifacts\n',...
		' 7. Split data\n',...
		' 9. Save and exit\n',...
		' 0. Exit without saving\n',...
		'Enter a menu item number: '];
	main_str = input(prompt,'s');
	if main_str == '0' 
		figure(gui); clf; break;

	elseif main_str == '9'
		figure(gui); clf; EEG = temp_EEG; EEG = eeg_checkset(EEG);	eeglab redraw; break;
		
	elseif main_str == '7'
		main_str = '';
		[temp_EEG, EEG, gui] = split(temp_EEG, EEG, gui);

	elseif main_str == '6'
		main_str = '';
		[temp_EEG, EEG, gui] = viz_arts(temp_EEG, EEG, gui);

	elseif main_str == '5'
		main_str = '';
		[temp_EEG, EEG, gui] = rm_chans_corr(temp_EEG, EEG, gui);

	elseif main_str == '4'
		main_str = '';
		[temp_EEG, EEG, gui] = quickplot(temp_EEG, EEG, gui);

	elseif main_str == '3'
		main_str = '';
		[temp_EEG, EEG, gui] = rm_amps(temp_EEG, EEG, gui);

	elseif main_str == '2'
		main_str = '';
		[temp_EEG, EEG, gui] = rm_sds(temp_EEG, EEG, gui);

	elseif main_str == '1'
		main_str = '';
		EEG = pop_loadset();
		temp_EEG = EEG;
		eeglab redraw;
	else
		disp(['Invalid input. You must type a number.' newline]);
	end
end