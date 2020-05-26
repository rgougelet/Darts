function eog = get_bipolar_eog(EEG)	

	uveog_i= find(strcmp({EEG.chanlocs.labels},'UVEOG'));
	lveog_i = find(strcmp({EEG.chanlocs.labels},'LVEOG'));
	lheog_i = find(strcmp({EEG.chanlocs.labels},'LHEOG'));
	rheog_i = find(strcmp({EEG.chanlocs.labels},'RHEOG'));
	
	uveog = EEG.data(uveog_i,:); %#ok<*FNDSB>
	lveog = EEG.data(lveog_i,:);
	lheog = EEG.data(lheog_i,:);
	rheog = EEG.data(rheog_i,:);
	veog = uveog-lveog;
	heog = lheog-rheog;
	eog = [veog;heog];
end