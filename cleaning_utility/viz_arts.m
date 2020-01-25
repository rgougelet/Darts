function [temp_EEG, EEG, gui] = viz_arts(temp_EEG, EEG, gui)
		figure(gui)
		vis_artifacts(temp_EEG,EEG, 'WindowLength', 40)
end

