%%
clc
load_r= [];
labs = {};
i = 0;
for subj = pipeline.Subjects
	for clust = subj.Clusters
		for epoch = clust.Epochs
			i = i+1;
			load_r(i).Subject = string(epoch.SubjectIndex);
			load_r(i).TargetAbsent = epoch.SNAP.isAbs;
			load_r(i).TargetPresent = epoch.SNAP.isPres;
			load_r(i).Distance = epoch.SNAP.Distance;
			load_r(i).DelayTime = epoch.DelayCondition;
			load_r(i).ThrowTime = epoch.SNAP.ThrowTime;
			load_r(i).FrontCluster = logical(epoch.ClusterIndex==1);
			load_r(i).BackCluster = logical(epoch.ClusterIndex==2);
			load_r(i).DelayAlphaAmp = epoch.BaselineCorrected.Delay.Alpha.Avg;
			load_r(i).DelayThetaAmp = epoch.BaselineCorrected.Delay.Theta.Avg;
			load_r(i).PreThrowAlphaAmp =epoch.BaselineCorrected.PreThrow.Alpha.Avg;
			load_r(i).PreThrowThetaAmp =epoch.BaselineCorrected.PreThrow.Theta.Avg;
			load_r(i).TDelayAlphaAmp = epoch.BaselineCorrected.Delay.Alpha.Avg*epoch.DelayCondition;
			load_r(i).TDelayThetaAmp = epoch.BaselineCorrected.Delay.Theta.Avg*epoch.DelayCondition;
			load_r(i).TPreThrowAlphaAmp =epoch.BaselineCorrected.PreThrow.Alpha.Avg*epoch.SNAP.ThrowTime;
			load_r(i).TPreThrowThetaAmp =epoch.BaselineCorrected.PreThrow.Theta.Avg*epoch.SNAP.ThrowTime;		
		end
	end
end
load_r = struct2table(load_r);
load_r=load_r(~any(ismissing(load_r),2),:);
rm_rows = abs(load_r.DelayThetaAmp)>15 & abs(load_r.DelayAlphaAmp)>15 & load_r.ThrowTime>4;
load_r(rm_rows,:) = [] ; % remove all nans
load_r.Distance = load_r.Distance*6.55; % convert to real-life distance

%% Tests for differences from baseline
% anterior
% target absent
r = avg_subj_r(load_r, 'FrontCluster','TargetAbsent');
boot_ttest(r.DelayThetaAmp,100000);
boot_ttest(r.DelayAlphaAmp,100000);
boot_ttest(r.PreThrowThetaAmp,100000);
boot_ttest(r.PreThrowAlphaAmp,100000);
% target present
r = avg_subj_r(load_r, 'FrontCluster','TargetPresent');
boot_ttest(r.DelayThetaAmp,100000);
boot_ttest(r.DelayAlphaAmp,100000);
boot_ttest(r.PreThrowThetaAmp,100000);
boot_ttest(r.PreThrowAlphaAmp,100000);
% time adjusted target absent
r = avg_subj_r(load_r, 'FrontCluster','TargetAbsent');
boot_ttest(r.TDelayThetaAmp,100000);
boot_ttest(r.TDelayAlphaAmp,100000);
boot_ttest(r.TPreThrowThetaAmp,100000);
boot_ttest(r.TPreThrowAlphaAmp,100000);
% time adjusted target present
r = avg_subj_r(load_r, 'FrontCluster','TargetPresent');
boot_ttest(r.TDelayThetaAmp,100000);
boot_ttest(r.TDelayAlphaAmp,100000);
boot_ttest(r.TPreThrowThetaAmp,100000);
boot_ttest(r.TPreThrowAlphaAmp,100000);
% posterior
% target absent
r = avg_subj_r(load_r, 'BackCluster','TargetAbsent');
boot_ttest(r.DelayThetaAmp,100000);
boot_ttest(r.DelayAlphaAmp,100000);
boot_ttest(r.PreThrowThetaAmp,100000);
boot_ttest(r.PreThrowAlphaAmp,100000);
% target present
r = avg_subj_r(load_r, 'BackCluster','TargetPresent');
boot_ttest(r.DelayThetaAmp,100000);
boot_ttest(r.DelayAlphaAmp,100000);
boot_ttest(r.PreThrowThetaAmp,100000);
boot_ttest(r.PreThrowAlphaAmp,100000);
% time adjusted target absent
r = avg_subj_r(load_r, 'BackCluster','TargetAbsent');
boot_ttest(r.TDelayThetaAmp,100000);
boot_ttest(r.TDelayAlphaAmp,100000);
boot_ttest(r.TPreThrowThetaAmp,100000);
boot_ttest(r.TPreThrowAlphaAmp,100000);
% time adjusted target present
r = avg_subj_r(load_r, 'BackCluster','TargetPresent');
boot_ttest(r.TDelayThetaAmp,100000);
boot_ttest(r.TDelayAlphaAmp,100000);
boot_ttest(r.TPreThrowThetaAmp,100000);
boot_ttest(r.TPreThrowAlphaAmp,100000);

%% Tests for differences between Target Absent and Target Present
% anterior
abs_r = avg_subj_r(load_r, 'FrontCluster','TargetAbsent');
pres_r = avg_subj_r(load_r, 'FrontCluster','TargetPresent');
perm_ttest(abs_r.DelayThetaAmp,pres_r.DelayThetaAmp,100000);
perm_ttest(abs_r.DelayAlphaAmp,pres_r.DelayAlphaAmp,100000);
perm_ttest(abs_r.PreThrowThetaAmp,pres_r.PreThrowThetaAmp,100000);
perm_ttest(abs_r.PreThrowAlphaAmp,pres_r.PreThrowAlphaAmp,100000);
% time adjusted
perm_ttest(abs_r.TDelayThetaAmp,pres_r.TDelayThetaAmp,100000);
perm_ttest(abs_r.TDelayAlphaAmp,pres_r.TDelayAlphaAmp,100000);
perm_ttest(abs_r.TPreThrowThetaAmp,pres_r.TPreThrowThetaAmp,100000);
perm_ttest(abs_r.TPreThrowAlphaAmp,pres_r.TPreThrowAlphaAmp,100000);
% posterior
abs_r = avg_subj_r(load_r, 'BackCluster','TargetAbsent');
pres_r = avg_subj_r(load_r, 'BackCluster','TargetPresent');
perm_ttest(abs_r.DelayThetaAmp,pres_r.DelayThetaAmp,100000);
perm_ttest(abs_r.DelayAlphaAmp,pres_r.DelayAlphaAmp,100000);
perm_ttest(abs_r.PreThrowThetaAmp,pres_r.PreThrowThetaAmp,100000);
perm_ttest(abs_r.PreThrowAlphaAmp,pres_r.PreThrowAlphaAmp,100000);
% time adjusted
perm_ttest(abs_r.TDelayThetaAmp,pres_r.TDelayThetaAmp,100000);
perm_ttest(abs_r.TDelayAlphaAmp,pres_r.TDelayAlphaAmp,100000);
perm_ttest(abs_r.TPreThrowThetaAmp,pres_r.TPreThrowThetaAmp,100000);
perm_ttest(abs_r.TPreThrowAlphaAmp,pres_r.TPreThrowAlphaAmp,100000);

%% Correlations with distance
% anterior
% target absent
r = avg_subj_r(load_r, 'FrontCluster','TargetAbsent');
perm_corr(r.Distance,r.DelayThetaAmp,100000);
perm_corr(r.Distance,r.DelayAlphaAmp,100000);
perm_corr(r.Distance,r.PreThrowThetaAmp,100000);
perm_corr(r.Distance,r.PreThrowAlphaAmp,100000);
% target present
r = avg_subj_r(load_r, 'FrontCluster','TargetPresent');
perm_corr(r.Distance,r.DelayThetaAmp,100000);
perm_corr(r.Distance,r.DelayAlphaAmp,100000);
perm_corr(r.Distance,r.PreThrowThetaAmp,100000);
perm_corr(r.Distance,r.PreThrowAlphaAmp,100000);
% time adjusted target absent
r = avg_subj_r(load_r, 'FrontCluster','TargetAbsent');
perm_corr(r.Distance,r.TDelayThetaAmp,100000);
perm_corr(r.Distance,r.TDelayAlphaAmp,100000);
perm_corr(r.Distance,r.TPreThrowThetaAmp,100000);
perm_corr(r.Distance,r.TPreThrowAlphaAmp,100000);
% time adjusted target present
r = avg_subj_r(load_r, 'FrontCluster','TargetPresent');
perm_corr(r.Distance,r.TDelayThetaAmp,100000);
perm_corr(r.Distance,r.TDelayAlphaAmp,100000);
perm_corr(r.Distance,r.TPreThrowThetaAmp,100000);
perm_corr(r.Distance,r.TPreThrowAlphaAmp,100000);

% posterior
% target absent
r = avg_subj_r(load_r, 'BackCluster','TargetAbsent');
perm_corr(r.Distance,r.DelayThetaAmp,100000);
perm_corr(r.Distance,r.DelayAlphaAmp,100000);
perm_corr(r.Distance,r.PreThrowThetaAmp,100000);
perm_corr(r.Distance,r.PreThrowAlphaAmp,100000);
% target present
r = avg_subj_r(load_r, 'BackCluster','TargetPresent');
perm_corr(r.Distance,r.DelayThetaAmp,100000);
perm_corr(r.Distance,r.DelayAlphaAmp,100000);
perm_corr(r.Distance,r.PreThrowThetaAmp,100000);
perm_corr(r.Distance,r.PreThrowAlphaAmp,100000);
% time adjusted target absent
r = avg_subj_r(load_r, 'BackCluster','TargetAbsent');
perm_corr(r.Distance,r.TDelayThetaAmp,100000);
perm_corr(r.Distance,r.TDelayAlphaAmp,100000);
perm_corr(r.Distance,r.TPreThrowThetaAmp,100000);
perm_corr(r.Distance,r.TPreThrowAlphaAmp,100000);
% time adjusted target present
r = avg_subj_r(load_r, 'BackCluster','TargetPresent');
perm_corr(r.Distance,r.TDelayThetaAmp,100000);
perm_corr(r.Distance,r.TDelayAlphaAmp,100000);
perm_corr(r.Distance,r.TPreThrowThetaAmp,100000);
perm_corr(r.Distance,r.TPreThrowAlphaAmp,100000);

%% Correlations with throw time
% anterior
% target absent
r = avg_subj_r(load_r, 'FrontCluster','TargetAbsent');
perm_corr(r.ThrowTime,r.DelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.DelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowAlphaAmp,100000);
% target present
r = avg_subj_r(load_r, 'FrontCluster','TargetPresent');
perm_corr(r.ThrowTime,r.DelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.DelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowAlphaAmp,100000);
% time adjusted target absent
r = avg_subj_r(load_r, 'FrontCluster','TargetAbsent');
perm_corr(r.ThrowTime,r.TDelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.TDelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowAlphaAmp,100000);
% time adjusted target present
r = avg_subj_r(load_r, 'FrontCluster','TargetPresent');
perm_corr(r.ThrowTime,r.TDelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.TDelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowAlphaAmp,100000);

% posterior
% target absent
r = avg_subj_r(load_r, 'BackCluster','TargetAbsent');
perm_corr(r.ThrowTime,r.DelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.DelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowAlphaAmp,100000);
% target present
r = avg_subj_r(load_r, 'BackCluster','TargetPresent');
perm_corr(r.ThrowTime,r.DelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.DelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.PreThrowAlphaAmp,100000);
% time adjusted target absent
r = avg_subj_r(load_r, 'BackCluster','TargetAbsent');
perm_corr(r.ThrowTime,r.TDelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.TDelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowAlphaAmp,100000);
% time adjusted target present
r = avg_subj_r(load_r, 'BackCluster','TargetPresent');
perm_corr(r.ThrowTime,r.TDelayThetaAmp,100000);
perm_corr(r.ThrowTime,r.TDelayAlphaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowThetaAmp,100000);
perm_corr(r.ThrowTime,r.TPreThrowAlphaAmp,100000);