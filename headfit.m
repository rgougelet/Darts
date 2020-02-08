function EEG = headfit(EEG,subjn)
	if subjn == '607'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[1.5 -17 -5.5 -0.066831 -0.12526 -1.6351 1.12 1.2 1.24] ,'chansel',[1:128] );
	end
	if subjn == '619'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.44856 -22.3762 -11.4453 -0.056364 0.045909 -1.5687 1.0518 1.2099 1.242] ,'chansel',[1:128] );
	end
	if subjn == '608'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.25 -17 -11 -0.052681 -0.027814 -1.5364 1.1 1.14 1.3294] ,'chansel',[1:128] );
	end
	if subjn == '616'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-1.1742 -24.3343 -5.6142 -0.093101 -0.076213 -1.5722 1.0524 1.2382 1.1889] ,'chansel',[1:128] );
	end
	if subjn == '627'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.65647 -18.4325 -4.012 -0.11676 -0.052117 -1.5699 1.1455 1.259 1.2161] ,'chansel',[1:128] );
	end
	if subjn == '580'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[1.75 -17 -2.5 -0.047721 -0.033688 -1.5893 1.075 1 1.14] ,'chansel',[1:128] );
	end
	if subjn == '579'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.5 -14.5 -18 -0.026591 0.0092933 -1.5633 1.12 1.135 1.375] ,'chansel',[1:128] );
	end
	if subjn == '571'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[2 -16 -9 -0.1872 0.064588 -1.5753 1.15 1.2 1.3] ,'chansel',[1:128] );
	end
	if subjn == '621'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.59856 -18.1779 -10.8712 -0.12325 -0.039926 -1.5485 1.0944 1.1884 1.3] ,'chansel',[1:128] );
	end
	if subjn == '631'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.50053 -22.9542 -13.4548 -0.13076 0.034242 -1.566 1.256 1.2096 1.3705] ,'chansel',[1:128] );
	end
end

