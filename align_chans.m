clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
% script_dir = 'G:/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
eeglab nogui;
data_dir = [script_dir,'data/'];
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
%%
for subj_i = 1:length(subjs_to_include)
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'_eeg_512_trim_ic.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	EEG = pop_select( EEG,'nochannel', {'M1','M2','UVEOG', 'LVEOG', 'LHEOG', 'RHEOG'});
	
	%%
	close all
		% 	ref_locs = 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\elec\standard_1020.elc';
		ref_locs = [script_dir, 'eeglab/plugins/dipfit3.3/standard_BEM/elec/standard_1020.elc'];
		r = readlocs(ref_locs);
		r = r(4:93);
	
% 	% 	ref_locs = 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\elec\standard_1005.elc';
% 	ref_locs = [script_dir, 'eeglab/plugins/dipfit3.3/standard_BEM/elec/standard_1005.elc'];
% 	r = readlocs(ref_locs);
% 	r = r(4:end-4);
	
	lab_is = true(1,length({r.labels}));
	for lab_i = 1:length({r.labels})
		if any(strfind(r(lab_i).labels,num2str(10))) || any(strfind(r(lab_i).labels,num2str(9))) || any(strfind(r(lab_i).labels,'Iz'))
			lab_is(lab_i) = false;
		end
	end
	r = r(lab_is);
	r3 = num2cell(ch2three(r),1);

	% plot manual coreg
	mEEG = headfit(EEG,subj_id);
	mc = ch2three(mEEG.chanlocs);
	transfmat = traditionaldipfit(mEEG.dipfit.coord_transform);
	tmc = transfmat*[ mc ones(size(mc,1),1) ]';
	tmc = (tmc(1:3,:))';
	mc3 = num2cell(tmc,1);
	figure;
	plot3(mc3{:},'*'); hold on;
	plot3(r3{:},'*'); hold on;
	xlim([-100,100]); ylim([-150, 150]); zlim([-30,110]);

	%%
	EEG =	pop_chanedit(EEG, 'nosedir','-Y');
	
%%
	% plot original locations
	c3 = num2cell(ch2three(EEG.chanlocs),1);
	figure;	
	plot3(c3{:},'*'); hold on;
	plot3(r3{:},'*'); hold on;
	xlim([-100,100]); ylim([-150, 150]); zlim([-30,110]);
	
	%% plot transformed locations using pcreg
	cl = pointCloud(ch2three(EEG.chanlocs));
	rl = pointCloud(ch2three(r));
	[tform, ~] = pcregrigid(cl,rl,'Extrapolate',true);
	tcl = pctransform(cl,tform);
	tc3 = num2cell(tcl.Location,1);
	figure;
	plot3(tc3{:},'*'); hold on;
	plot3(r3{:},'*'); hold on;
	xlim([-100,100]); ylim([-150, 150]); zlim([-30,110]);

	%% plot transformed locations using pcreg, starting with manual
	mEEG =	pop_chanedit(mEEG, 'nosedir','-Y');
	cl = pointCloud(tmc);
	rl = pointCloud(ch2three(r));
	[tform, ~] = pcregrigid(cl,rl,'Extrapolate',true);
	tcl = pctransform(cl,tform);
	tc3 = num2cell(tcl.Location,1);
	figure;
	plot3(tc3{:},'*'); hold on;
	plot3(r3{:},'*'); hold on;
	xlim([-100,100]); ylim([-150, 150]); zlim([-30,110]);

end

%% trash
% 	EEG.chanlocs = three2ch(c,tcl.Location);
% 		[~, ~] = coregister(EEG.chanlocs,ref_locs,in_cell{:});

% 	tpc_locs;
% 	new_locs = EEG.chanlocs(1:128);
% 	x = num2cell(tpc_locs.Location(:,1));
% 	[new_locs.X] = x{:};
% 	y = num2cell(tpc_locs.Location(:,2));
% 	[new_locs.Y] = y{:};
% 	z = num2cell(tpc_locs.Location(:,3));
% 	[new_locs.Z] = z{:};
% 	new_locs.Y = tpc_locs.Location(:,2);
% 	new_locs.Z = tpc_locs.Location(:,3);

% 	TMP          = eeg_emptyset;
% 	TMP.chanlocs = readlocs(EEG.chanlocs);
% 	TMP.chaninfo = EEG.chaninfo;
% 	TMP.nbchan = length(TMP.chanlocs);
% 	cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
% 	elec1 = cfg.elec;
% 	if isfield(elec1, 'elecpos')
% 		elec1.pnt = elec1.elecpos;
% 	end
% 	
% 	TMP          = eeg_emptyset;
% 	TMP.chanlocs = readlocs(ref_locs);
% 	TMP.nbchan = length(TMP.chanlocs);
% 	cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
% 	elec2 = cfg.elec;
% 	if isfield(elec2, 'elecpos')
% 		elec2.pnt = elec2.elecpos;
% 	end
	
% 	refs = {...
% 		'F1'
% 		'Fz'
% 		'F2'
% 		'FC1'
% 		'FCz'
% 		'FC2'
% 		'C1'
% 		'Cz'
% 		'C2'
% 		'CP1'
% 		'CPz'
% 		'CP2'
% 		'P1'
% 		'Pz'
% 		'P2'
% 		};
% 	[~,~,ref_i] = intersect(refs, elec2.label);
% 		elec2.pnt = (elec2.pnt(ref_i,:));
	
% 	elec1.pnt = (elec1.pnt(1:128,:));
% 	pnt1 = elec1.pnt;
% 	elec2.label = (elec2.label(4:93));
% 	elec2.pnt = (elec2.pnt(4:93,:));
% 	lab_is = true(1,numel(elec2.label));
% 	for lab_i = 1:length(elec2.label)
% 		if any(strfind(elec2.label{lab_i},num2str(10))) || any(strfind(elec2.label{lab_i},num2str(9))) || any(strfind(elec2.label{lab_i},'Iz'))
% 			lab_is(lab_i) = false;
% 		end
% 	end
% 	elec2.label = (elec2.label(lab_is));
% 	elec2.pnt = (elec2.pnt(lab_is,:));
% 	pnt2 = (elec2.pnt);
% 	init_t = EEG.dipfit.coord_transform;
	
% 	transfmat = traditionaldipfit(init_t);
% 	pnt1 = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
% 	pnt1 = (pnt1(1:3,:))';
% 	pc_elec = pointCloud(pnt1);
% 	pc_ref = pointCloud(pnt2);
% 	pcshow(ptCloud_elec)
% 	tform = pcregrigid(pc_elec,pc_ref,'Extrapolate',true);
% 	fit_t = tform.T;
% 	chans = EEG.chanlocs(1:128);
% 	locs = [chans.X;chans.Y;chans.Z]';
% 	pc_locs = pointCloud(locs);
% 	tpc_locs = pctransform(pc_locs,tform);
% 	tpc_locs;
% 	new_locs = EEG.chanlocs(1:128);
% 	x = num2cell(tpc_locs.Location(:,1));
% 	[new_locs.X] = x{:};
% 	y = num2cell(tpc_locs.Location(:,2));
% 	[new_locs.Y] = y{:};
% 	z = num2cell(tpc_locs.Location(:,3));
% 	[new_locs.Z] = z{:};
% 	new_locs.Y = tpc_locs.Location(:,2);
% 	new_locs.Z = tpc_locs.Location(:,3);

% 	t1s = -4:1:4;
% 	t2s = -30:1:0;
% 	t3s = 10:1:30;
% 	t4s = 1:0.05:1.5;
% 	t5s = 1:0.05:1.5;
% 	t6s = 1:0.05:1.5;
% 	t1s = -4:1:4;
% 	t2s = -30:1:0;
% 	t3s = -20:1:20;
% 	t4s = init_t(7);
% 	t5s = .9:0.025:1.6;
% 	t6s = .9:0.025:1.6;
% 	numts = length(t1s)*length(t2s)*length(t3s)*length(t4s)*length(t5s)*length(t6s);
% 	
% 	ts = nan(numts,9);
% 	t_i = 0;
% 	for t1 = t1s
% 		for t2 = t2s
% 			for t3 = t3s
% 				for t4 = t4s
% 					for t5 = t5s
% 						for t6 = t6s
% 							t_i = t_i+1;
% 							ts(t_i,:) = [t1 t2 t3 init_t(4) init_t(5) init_t(6) t4 t5 t6];
% 						end
% 					end
% 				end
% 			end
% 		end
% 	end
	
	% get distances for all transformations
% 	ds = nan(10,numts);
% 	parfor d_i = 1:numts
% 		t = ts(d_i,:);
% 		transfmat = traditionaldipfit(t);
% 		pnt1 = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
% 		pnt1 = (pnt1(1:3,:))';
% 		
% 		r1 = repmat(pnt1,length(pnt2),1);
% 		r2 = repelem(pnt2,length(pnt1),1);
% 		r = sum((r1-r2).^2,2);
% 		
% 		d = median(r);
% 		
% 		ds(:,d_i) = [t, d];
% 	end
	
% 	[min_d,i] = min(ds(10,:));
% 	min_t = ds(1:9,i)';
% 	
% 	transfmat = traditionaldipfit(min_t);
% 	pnt1 = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
% 	pnt1 = (pnt1(1:3,:))';
% 	r1 = repmat(pnt1,length(pnt2),1);
% 	r2 = repelem(pnt2,length(pnt1),1);
% 	r = sum((r1-r2).^2,2);
% 	init_d = median(r)
% 	
% 	transfmat = traditionaldipfit(init_t);
% 	pnt1 = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
% 	pnt1 = (pnt1(1:3,:))';
% 	
% 	r1 = repmat(pnt1,length(pnt2),1);
% 	r2 = repelem(pnt2,length(pnt1),1);
% 	r = sum((r1-r2).^2,2);
% 	min_d = median(r)

% plot original manual transformation
% 	in_cell = {'mesh',...
% 		'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_vol.mat',...
% 		'transform',...
% 		[],...
% 		'chaninfo1',...
% 		EEG.chaninfo,...
% 		'helpmsg',...
% 		'on',...
% 		'transform',...
% 		init_t
% 		};
% 	[~, ~] = coregister(EEG.chanlocs,ref_locs,in_cell{:});
	
	% plot brute-force matched transformation
% 	in_cell = {'mesh',...
% 		'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_vol.mat',...
% 		'transform',...
% 		[],...
% 		'chaninfo1',...
% 		EEG.chaninfo,...
% 		'helpmsg',...
% 		'on',...
% 		'transform',...
% 		init_t*fit_t
% 		};
% EEG.chaninfo.nosedir = '-Y';
% 	in_cell = {'mesh',...
% 		'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_vol.mat',...
% 		'transform',...
% 		[],...
% 		'chaninfo1',...
% 		EEG.chaninfo,...
% 		'helpmsg',...
% 		'on'
% 		};
% 	
% 	[~, ~] = coregister(new_locs,ref_locs,in_cell{:});