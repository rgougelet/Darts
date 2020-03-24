for c_i = 1:length(EEG.dipfit.model)
	pos1 = EEG.dipfit.model(c_i).posxyz(1,:);
	pos2 = EEG.dipfit.model(c_i).posxyz(2,:);
	pos1s(c_i,:) = pos1;
	pos2s(c_i,:) = pos2;

end
%%

% plot3(pos1s(:,1),pos1s(:,2),pos1s(:,3), 'r*')
%%
% plot3(pos2s(:,1),pos2s(:,2),pos2s(:,3), 'b*')

dists = sum(pos1s.^2,2);

% figure; hist(dists,100)
[s, i] = sort(dists,'descend')
%%
EEG = pop_iclabel(EEG, default);
lab = EEG.etc.ic_classification.ICLabel;
[~, sort_i] = sort(lab.classifications(:,1),'ascend')
eegplot(EEG.icaact(sort_i,:)./std(EEG.icaact(sort_i,:),[],2), 'dispchans',32, 'limits', [0 5])
% eegplot(EEG.icaact(sort_i,:), 'dispchans',32, 'limits', [0 5])

%%
[sources X Y Z XE YE ZE] = dipplot(EEG.dipfit.model(i),...
 'meshdata', 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_vol.mat',...
 'coordformat','MNI',...
 'mri', 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_mri.mat',...
 'normlen', 'on');

%%
sources = [];
sources(1).posxyz = [0 0 0];   % position for the first dipole
 sources(1).momxyz = [  0 58 -69];   % orientation for the first dipole
 sources(1).rv     = 0.036; 

[sources X Y Z XE YE ZE] = dipplot(sources,...
 'meshdata', 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_vol.mat',...
 'coordformat','MNI',...
 'mri', 'G:\darts\eeglab\plugins\dipfit3.3\standard_BEM\standard_mri.mat',...
 'normlen', 'on',...
 'projlines', 'on');