clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab/')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'deps/'])
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

for subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir, subj_id,'*_lab.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	old_EEG = EEG;

%% find correlated frequencies, i.e. freq-freq coupling
	inst_freqs = instfreq(EEG.icaact(:,:)',EEG.srate,'method', 'hilbert');
	%%
	ff_corr = corr(inst_freqs);
	ff_corr(tril(true(size(ff_corr)))) = nan;

	for c1 = 1:134
		[m, c2] = max(ff_corr(c1,:));
		ms(c1) = m;
		cs(c1,:) = [c1,c2];
	end
	[~, sort_i] = sort(ms,'ascend');
	cs = cs(sort_i,:);
% 	d1 = EEG.icaact(cs(:,1),:);
	d1 = inst_freqs(:,cs(:,1))';
	d1 = d1./std(d1,[],2);
% 	d2 = EEG.icaact(cs(:,2),:);
	d2 = inst_freqs(:,cs(:,2))';
	d2 = d2./std(d2,[],2);
	eegplot(d1, 'data2', d2, 'dispchans',32, 'winlength',5);

% 	[~,m1] = max(ff_corr,[],1)
% 	[~,m2] = max(ff_corr,[],2)
% 	maxs = [m1;m2']';
% 	sort_i = sort(ff_corr(meshgrid(m1,m2)));
% 	inds = nchoosek(1:134,2);
% 	lin_inds = sub2ind(size(ff_corr),inds(:,1),inds(:,2));
% 	lin_ff_cor = ff_corr(lin_inds);
% 	[~,sort_i] = sort(ff_corr(lin_inds),'descend');
% 	[sort_1, sort_2] = ind2sub(size(ff_corr),lin_inds(sort_i));
% 	for i = 1:length(sort_i)
% 		corrs(i) = ff_corr(sort_1(i),sort_2(i));
% 	end
% 	eegplot(EEG.icaact(m1(sort_i),:), 'data2', EEG.icaact(m2(sort_i),:));

% 	ff_corr(tril(true(size(ff_corr)))) = nan;
% 	imagesc(ff_corr);
% 	[~, I] = nanmax(ff_corr);
% 	[i,j] = sub2ind(size(ff_corr),I);
% 	eegplot(EEG.icaact(i,:), 'data2', EEG.icaact(j,:));
% 	n = size(EEG.icaact,1);
% 	max(ff_corr)
% 	kurts = kurtosis(inst_freqs);
% 	[sort_c, sort_i] = sort(kurts,'ascend');
% 	eegplot(EEG.icaact(sort_i,:), 'dispchans',32, 'limits', [0 5], 'winlength',5)
% 	eegplot(inst_freqs', 'dispchans',32, 'limits', [0 5], 'winlength',5)
%% find ICs whose frequency oscillates
% 	inst_df = instfreq(inst_freqs,EEG.srate,'method', 'hilbert');
% 	kurts = kurtosis(inst_ff);
% 	[sort_c, sort_i] = sort(kurts,'ascend');
% 	eegplot(EEG.icaact(sort_i,:), 'dispchans',32, 'winlength',5)
% 	eegplot(inst_ff(:,sort_i)', 'dispchans',32, 'winlength',5)

end