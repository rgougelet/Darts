clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
eeglab nogui;
data_dir = [script_dir,'data/32/'];
addpath(data_dir)
subjs_to_include = {
	'571'
	'579'
	'580'
	'607'
	'608'
	'616'
	'619'
	'621'
	'627'
	'631'
	};
new_srate = 512;

for subj_i = 1%:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'_eeg_32.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	EEG = headfit(EEG,subj_id);
	ref_locs = 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\elec\standard_1020.elc';
	%%
	TMP          = eeg_emptyset;
	TMP.chanlocs = readlocs(EEG.chanlocs);
	TMP.chaninfo = EEG.chaninfo;
	TMP.nbchan = length(TMP.chanlocs);
	cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
	elec1 = cfg.elec;
	if isfield(elec1, 'elecpos')
		elec1.pnt = elec1.elecpos;
	end
	
	TMP          = eeg_emptyset;
	TMP.chanlocs = readlocs(ref_locs);
	TMP.nbchan = length(TMP.chanlocs);
	cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
	elec2 = cfg.elec;
	if isfield(elec2, 'elecpos')
		elec2.pnt = elec2.elecpos;
	end
	
	elec1.pnt = single(elec1.pnt(1:128,:));
	elec2.pnt = single(elec2.pnt(4:end-4,:));
	pnt2 = (elec2.pnt);
	
	t1s = -4:1:4;
	t2s = -30:1:0;
	t3s = 10:-1:-20;
	t4s = 1:0.05:1.3;
	t5s = 1:0.05:1.3;
	t6s = 1:0.05:1.4;
	numts = length(t1s)*length(t2s)*length(t3s)*length(t4s)*length(t5s)*length(t6s);
	init_t = EEG.dipfit.coord_transform;

	ts = nan(numts,9);
	t_i = 0;
	for t1 = t1s
		for t2 = t2s
			for t3 = t3s
				for t4 = t4s
					for t5 = t5s
						for t6 = t6s
							t_i = t_i+1;
							ts(t_i,:) = [t1 t2 t3 init_t(4) init_t(5) init_t(6) t4 t5 t6];
						end
					end
				end
			end
		end
	end

	ds = nan(10,numts);
	parfor d_i = 1:numts
		t = ts(d_i,:);
		transfmat = traditionaldipfit(t);
		pnt1 = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
		pnt1 = (pnt1(1:3,:))';

		r1 = repmat(pnt1,length(pnt2),1);
		r2 = repelem(pnt2,length(pnt1),1);
		d = sum(sum((r1-r2).^2));

		ds(:,d_i) = [t, d];
	end

	[~,i] = min(ds(10,:));
	min_t = ds(1:9,i)';
	
	in_cell = {'mesh',...
		'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_vol.mat',...
		'transform',...
		[],...
		'chaninfo1',...
		EEG.chaninfo,...
		'helpmsg',...
		'on',...
		'transform',...
		min_t
		};
	[elec1, ~] = coregister(EEG.chanlocs,ref_locs,in_cell{:});
	
end
% new_EEG = pop_dipfit_settings( EEG,...
%  'hdmfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_vol.mat',...
% 'coordformat','MNI',...
% 'mrifile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\standard_mri.mat',...
% 'chanfile','G:\\darts\\eeglab\\plugins\\dipfit3.3\\standard_BEM\\elec\\standard_1005.elc',...
% 'coord_transform',transform ,...
% 'chansel',[1:128] );
% x = [[1.5 -17 -5.5 -0.066831 -0.12526 -1.6351 1.12 1.2 1.24];
% [-0.44856 -22.3762 -11.4453 -0.056364 0.045909 -1.5687 1.0518 1.2099 1.242];
% [0.25 -17 -11 -0.052681 -0.027814 -1.5364 1.1 1.14 1.3294] ;
% [-1.1742 -24.3343 -5.6142 -0.093101 -0.076213 -1.5722 1.0524 1.2382 1.1889] ;
% [-0.65647 -18.4325 -4.012 -0.11676 -0.052117 -1.5699 1.1455 1.259 1.2161] ;
% [1.75 -17 -2.5 -0.047721 -0.033688 -1.5893 1.075 1 1.14] ;
% [0.5 -14.5 -18 -0.026591 0.0092933 -1.5633 1.12 1.135 1.375] ;
% [2 -16 -9 -0.1872 0.064588 -1.5753 1.15 1.2 1.3] ;
% [0.59856 -18.1779 -10.8712 -0.12325 -0.039926 -1.5485 1.0944 1.1884 1.3] ;
% [-0.50053 -22.9542 -13.4548 -0.13076 0.034242 -1.566 1.256 1.2096 1.3705] ];