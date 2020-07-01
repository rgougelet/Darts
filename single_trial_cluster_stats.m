function cluster_stats(pipeline,cluster_i)
%% Get single trial data
	i = 0;
	a_i = 0;
	pres_i = 0;
	for subj = pipeline.Subjects
		for clust = subj.Clusters(cluster_i)
			for epoch = clust.Epochs
				if ~strcmp(epoch.TrialString,epoch.SNAP.String); error('String mismatch'); end
				i = i+1;
				if epoch.SNAP.isAbs
					a_i = a_i+1;
					a.Delay.Theta(a_i,:) = epoch.BaselineCorrected.Delay.Theta.Avg;
					a.Delay.Alpha(a_i,:) = epoch.BaselineCorrected.Delay.Alpha.Avg;
					a.PreThrow.Theta(a_i,:) = epoch.BaselineCorrected.PreThrow.Theta.Avg;
					a.PreThrow.Alpha(a_i,:) = epoch.BaselineCorrected.PreThrow.Alpha.Avg;
					a.Delay.Cond(a_i,:) = epoch.SNAP.Delay; % snap timing seems better
					a.Distance(a_i,:) = epoch.SNAP.Distance;
					a.Throw.Time(a_i,:) = epoch.SNAP.ThrowTime;
				end
				if epoch.SNAP.isPres
					pres_i = pres_i+1;
					pres.Delay.Theta(pres_i,:) = epoch.BaselineCorrected.Delay.Theta.Avg;
					pres.Delay.Alpha(pres_i,:) = epoch.BaselineCorrected.Delay.Alpha.Avg;
					pres.PreThrow.Theta(pres_i,:) = epoch.BaselineCorrected.PreThrow.Theta.Avg;
					pres.PreThrow.Alpha(pres_i,:) = epoch.BaselineCorrected.PreThrow.Alpha.Avg;
					pres.Delay.Cond(pres_i,:) = epoch.SNAP.Delay;
					pres.Distance(pres_i,:) = epoch.SNAP.Distance;
					pres.Throw.Time(pres_i,:) = epoch.SNAP.ThrowTime;
				end
			end
		end
	end
clust.Label

%% Tests for differences from baseline
a.rmo = a.Throw.Time < 4;
pres.rmo = pres.Throw.Time < 4;
% target absent
boot_ttest(a.Delay.Theta,100000);
boot_ttest(a.Delay.Alpha,100000);
boot_ttest(a.PreThrow.Theta,100000);
boot_ttest(a.PreThrow.Alpha,100000);
% target present
boot_ttest(pres.Delay.Theta,100000);
boot_ttest(pres.Delay.Alpha,100000);
boot_ttest(pres.PreThrow.Theta,100000);
boot_ttest(pres.PreThrow.Alpha,100000);

% time adjusted target absent
boot_ttest(a.Delay.Theta.*a.Delay.Cond,100000);
boot_ttest(a.Delay.Alpha.*a.Delay.Cond,100000);
boot_ttest(a.PreThrow.Theta(a.rmo).*a.Throw.Time(a.rmo),100000);
boot_ttest(a.PreThrow.Alpha(a.rmo).*a.Throw.Time(a.rmo),100000);
% time adjusted target present
boot_ttest(pres.Delay.Theta.*pres.Delay.Cond,100000);
boot_ttest(pres.Delay.Alpha.*pres.Delay.Cond,100000);
boot_ttest(pres.PreThrow.Theta(pres.rmo).*pres.Throw.Time(pres.rmo),100000);
boot_ttest(pres.PreThrow.Alpha(pres.rmo).*pres.Throw.Time(pres.rmo),100000);

%% Tests for differences between Target Absent and Target Present
a.rmo = a.Throw.Time < 4;
pres.rmo = pres.Throw.Time < 4;

perm_ttest2(a.Delay.Theta,pres.Delay.Theta,'nperm',10000);
perm_ttest2(a.Delay.Alpha,pres.Delay.Alpha,'nperm',10000);
perm_ttest2(a.PreThrow.Theta,pres.PreThrow.Theta,'nperm',10000);
perm_ttest2(a.PreThrow.Alpha,pres.PreThrow.Alpha,'nperm',10000);
% time adjusted
perm_ttest2(a.Delay.Theta.*a.Delay.Cond,pres.Delay.Theta.*pres.Delay.Cond,'nperm',10000);
perm_ttest2(a.Delay.Alpha.*a.Delay.Cond,pres.Delay.Alpha.*pres.Delay.Cond,'nperm',10000);
perm_ttest2(a.PreThrow.Theta(a.rmo).*a.Throw.Time(a.rmo),...
pres.PreThrow.Theta(pres.rmo).*pres.Throw.Time(pres.rmo),'nperm',10000);
perm_ttest2(a.PreThrow.Alpha(a.rmo).*a.Throw.Time(a.rmo),...
pres.PreThrow.Alpha(pres.rmo).*pres.Throw.Time(pres.rmo),'nperm',10000);

%% Correlations with distance
a.rmo = a.Throw.Time < 4;
pres.rmo = pres.Throw.Time < 4;
% target absent
perm_corr(a.Distance, a.Delay.Theta); 
perm_corr(a.Distance, a.Delay.Alpha); 
perm_corr(a.Distance, a.PreThrow.Theta); 
perm_corr(a.Distance, a.PreThrow.Alpha); 
% target present
perm_corr(pres.Distance, pres.Delay.Theta); 
perm_corr(pres.Distance, pres.Delay.Alpha); 
perm_corr(pres.Distance, pres.PreThrow.Theta); 
perm_corr(pres.Distance, pres.PreThrow.Alpha); 

% time adjusted target absent
perm_corr(a.Distance, a.Delay.Theta.*a.Delay.Cond); 
perm_corr(a.Distance, a.Delay.Alpha.*a.Delay.Cond); 
perm_corr(a.Distance(a.rmo), a.PreThrow.Theta(a.rmo).*a.Throw.Time(a.rmo)); 
perm_corr(a.Distance(a.rmo), a.PreThrow.Alpha(a.rmo).*a.Throw.Time(a.rmo)); 
% time adjusted target present
perm_corr(pres.Distance, pres.Delay.Theta.*pres.Delay.Cond); 
perm_corr(pres.Distance, pres.Delay.Alpha.*pres.Delay.Cond); 
perm_corr(pres.Distance(pres.rmo), pres.PreThrow.Theta(pres.rmo).*pres.Throw.Time(pres.rmo)); 
perm_corr(pres.Distance(pres.rmo), pres.PreThrow.Alpha(pres.rmo).*pres.Throw.Time(pres.rmo)); 

%% Correlations with throw time
a.rmo = a.Throw.Time < 4;
pres.rmo = pres.Throw.Time < 4;
% target absent
perm_corr(a.Throw.Time(a.rmo), a.Delay.Theta(a.rmo));
perm_corr(a.Throw.Time(a.rmo), a.Delay.Alpha(a.rmo));
perm_corr(a.Throw.Time(a.rmo), a.PreThrow.Theta(a.rmo));
perm_corr(a.Throw.Time(a.rmo), a.PreThrow.Alpha(a.rmo));
% target present
perm_corr(pres.Throw.Time(pres.rmo), pres.Delay.Theta(pres.rmo));
perm_corr(pres.Throw.Time(pres.rmo), pres.Delay.Alpha(pres.rmo));
perm_corr(pres.Throw.Time(pres.rmo), pres.PreThrow.Theta(pres.rmo),100000);
perm_corr(pres.Throw.Time(pres.rmo), pres.PreThrow.Alpha(pres.rmo));

% target absent
perm_corr(a.Throw.Time(a.rmo), a.Delay.Theta(a.rmo).*a.Delay.Cond(a.rmo)); 
perm_corr(a.Throw.Time(a.rmo), a.Delay.Alpha(a.rmo).*a.Delay.Cond(a.rmo)); 

% target present
perm_corr(pres.Throw.Time(pres.rmo), pres.Delay.Theta(pres.rmo).*pres.Delay.Cond(pres.rmo)); 
perm_corr(pres.Throw.Time(pres.rmo), pres.Delay.Alpha(pres.rmo).*pres.Delay.Cond(pres.rmo)); 

%% Crude coupling which turns out to be correlation
a.rmo = a.Throw.Time < 4;
pres.rmo = pres.Throw.Time < 4;
% target absent coupling
boot_ttest(a.Delay.Theta.*a.Delay.Alpha,10000)
boot_ttest(a.PreThrow.Theta.*a.PreThrow.Alpha,10000)
% target present coupling
boot_ttest(pres.Delay.Theta.*pres.Delay.Alpha,10000)
boot_ttest(pres.PreThrow.Theta.*pres.PreThrow.Alpha,10000)

% time adjusted target absent coupling
boot_ttest(a.Delay.Theta.*a.Delay.Alpha.*a.Delay.Cond,10000)
boot_ttest(a.PreThrow.Theta.*a.PreThrow.Alpha.*a.Delay.Cond,10000)
% time adjusted target present coupling 
boot_ttest(pres.Delay.Theta.*pres.Delay.Alpha.*pres.Delay.Cond,10000)
boot_ttest(pres.PreThrow.Theta.*pres.PreThrow.Alpha.*pres.Delay.Cond,10000)

% coupling
perm_corr(a.Distance, a.Delay.Alpha.*a.Delay.Theta); 
perm_corr(a.Distance, a.PreThrow.Alpha.*a.PreThrow.Theta); 
perm_corr(a.Distance, a.Delay.Alpha.*a.Delay.Theta.*a.Delay.Cond); 
perm_corr(a.Distance(a.rmo), a.PreThrow.Alpha(a.rmo).*a.PreThrow.Theta(a.rmo).*a.Throw.Time(a.rmo)); 
perm_corr(pres.Distance, pres.Delay.Alpha.*pres.Delay.Theta); 
perm_corr(pres.Distance, pres.PreThrow.Alpha.*pres.PreThrow.Theta); 
perm_corr(pres.Distance, pres.Delay.Alpha.*pres.Delay.Theta.*pres.Delay.Cond); 
perm_corr(pres.Distance(pres.rmo), pres.PreThrow.Alpha(pres.rmo).*pres.PreThrow.Theta(pres.rmo).*pres.Throw.Time(pres.rmo)); 
