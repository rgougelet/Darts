function EEG = headfit(EEG,subjn)
	if subjn == '607'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0 -17 -4 -0.066831 -0.12526 -1.6351 1.12 1.16 1.23] ,'chansel',[1:128] );
	end
	if subjn == '619'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.5 -17 -12 -0.056364 0.045909 -1.5687 1.05 1.08 1.25] ,'chansel',[1:128] );
	end
	if subjn == '608'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.5 -16 -12 -0.052681 -0.027814 -1.5364 1.075 1.16 1.3294] ,'chansel',[1:128] );
	end
	if subjn == '616'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.5 -16 -5 -0.093101 -0.076213 -1.5722 1.025 1.14 1.17] ,'chansel',[1:128] );
	end
	if subjn == '627'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[1 -15.5 -6 -0.11676 -0.052117 -1.5699 1.1 1.35 1.25] ,'chansel',[1:128] );
	end
	if subjn == '580'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[1 -17 0.5 -0.047721 -0.033688 -1.5893 1.075 .99 1.09] ,'chansel',[1:128] );
	end
	if subjn == '579'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0 -15.5 -15 -0.026591 0.0092933 -1.5633 1.1075 1.075 1.34] ,'chansel',[1:128] );
	end
	if subjn == '571'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0 -14.5 -4.5 -0.1872 0.064588 -1.5753 1.1 1.13 1.25] ,'chansel',[1:128] );
	end
	if subjn == '621'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.59856 -16 -8 -0.12325 -0.039926 -1.5485 1.15 1.075 1.25] ,'chansel',[1:128] );
	end
	if subjn == '631'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.5 -15 -14 -0.13076 0.034242 -1.566 1.225 1.08 1.3705] ,'chansel',[1:128] );
	end
end

