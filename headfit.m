function EEG = headfit(EEG,subjn)
	if subjn == '607'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-1.7107 -19.6782 -1.8377 -0.066779 -0.12592 -1.6348 1.1396 1.3765 1.2028] ,'chansel',[1:128] );
	end
	if subjn == '619'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.44856 -22.3762 -11.4453 -0.056364 0.045909 -1.5687 1.0518 1.2099 1.242] ,'chansel',[1:128] );
	end
	if subjn == '608'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-2.0149 -20.2807 -12.6854 -0.052662 -0.027804 -1.5365 1.1009 1.294 1.3294] ,'chansel',[1:128] );
	end
	if subjn == '616'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-1.1742 -24.3343 -5.6142 -0.093101 -0.076213 -1.5722 1.0524 1.2382 1.1889] ,'chansel',[1:128] );
	end
	if subjn == '627'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.65647 -18.4325 -4.012 -0.11676 -0.052117 -1.5699 1.1455 1.259 1.2161] ,'chansel',[1:128] );
	end
	if subjn == '580'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[1.0542 -20.2478 -3.7656 -0.04772 -0.033688 -1.5893 1.1508 1.1508 1.1508] ,'chansel',[1:128] );
	end
	if subjn == '579'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','/data/common/matlab/eeglab/plugins/dipfit2.3/standard_BEM/standard_vol.mat','coordformat','MNI','mrifile','/data/common/matlab/eeglab/plugins/dipfit2.3/standard_BEM/standard_mri.mat','chanfile','/data/common/matlab/eeglab/plugins/dipfit2.3/standard_BEM/elec/standard_1005.elc','coord_transform',[6 -13 -37 -0.05 0 -1.575 1.1 1.12 1.18] ,'chansel', [1:128] );
	end
	if subjn == '571'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-2.8503 -20.3543 -6.0605 -0.14437 0.087909 -1.5728 1.1503 1.3409 1.3264] ,'chansel',[1:128] );	
	end
	if subjn == '621'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.59856 -18.1779 -10.8712 -0.12325 -0.039926 -1.5485 1.0944 1.1884 1.3] ,'chansel',[1:128] );
	end
	if subjn == '631'
		EEG = pop_dipfit_settings( EEG, 'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat','chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[-0.50053 -22.9542 -13.4548 -0.13076 0.034242 -1.566 1.256 1.2096 1.3705] ,'chansel',[1:128] );
	end
end

