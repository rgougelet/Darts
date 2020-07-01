function pipeline = get_behav_data(pipeline)
tic
% load XLSX SNAP data
[num,txt,~] = xlsread([pipeline.DataDir,'behavioral_data_reduced.xlsx']);
r = struct;
headers = txt(1,:);
for k=1:numel(headers)
	for ri = 1:length(num)
		r(ri).(headers{k})=num(ri,k);
	end
end
for subj_i = 1:pipeline.nSubjects
	pipeline.Subjects(subj_i).nTrials = sum([r.subject] == str2double(pipeline.Subjects(subj_i).ID));
	pipeline.Subjects(subj_i).SNAP.Indices = []; % remove at next run through
	pipeline.Subjects(subj_i).SNAP.Indices.Logical =	[r.subject] == str2double(pipeline.Subjects(subj_i).ID);
	pipeline.Subjects(subj_i).SNAP.Indices.Sequential =	find(pipeline.Subjects(subj_i).SNAP.Indices.Logical==1);

	subj = pipeline.Subjects(subj_i);
	for clust_i = 1:pipeline.nClusters
		clust = subj.Clusters(clust_i);
		snap_trial_strs = [...
			num2str([r(subj.SNAP.Indices.Logical).delay]'),...
			num2str([r(subj.SNAP.Indices.Logical).position]','%02d'),...
			num2str([r(subj.SNAP.Indices.Logical).pres]');];
		c = 1;
		for epoch = clust.Epochs
			if strcmp(subj.ID,'616') || strcmp(subj.ID,'621') || strcmp(subj.ID,'627')
				for snap_str_i = c:length(snap_trial_strs)
					if strcmp(epoch.TrialString,snap_trial_strs(snap_str_i,:))
						clust.Epochs(epoch.Index).TrialIndex = snap_str_i;
						clust.Epochs(epoch.Index).SNAP.Index = pipeline.Subjects(subj_i).SNAP.Indices.Sequential(snap_str_i);
						c = c + 1;
						break
					end
				end
			elseif strcmp(subj.ID, '580') % subj 580 missing first 10 trials
				clust.Epochs(epoch.Index).TrialIndex = epoch.Index+10;
				clust.Epochs(epoch.Index).SNAP.Index = pipeline.Subjects(subj_i).SNAP.Indices.Sequential(epoch.Index+10);
			else
				clust.Epochs(epoch.Index).TrialIndex = epoch.Index;
				clust.Epochs(epoch.Index).SNAP.Index = pipeline.Subjects(subj_i).SNAP.Indices.Sequential(epoch.Index);
			end
		end
		for epoch = clust.Epochs
			snap_i = epoch.SNAP.Index;
			epoch.SNAP.String = [...
			num2str([r(snap_i).delay]'),...
			num2str([r(snap_i).position]','%02d'),...
			num2str([r(snap_i).pres]')];
			epoch.SNAP.ThrowTime = r(snap_i).throwtime;
			epoch.SNAP.Delay = r(snap_i).delay;
			epoch.SNAP.Distance = r(snap_i).distance;
			epoch.SNAP.DelayTime = r(snap_i).delaystilltime;
			isPres= logical(r(snap_i).pres);
			epoch.SNAP.isPres = isPres;
			epoch.SNAP.isAbs = ~isPres;
			if isPres == 0; cond = 'Target Absent';
			elseif isPres == 1; cond = 'Target Present'; end
			epoch.SNAP.Condition = cond;
			clust.Epochs(epoch.Index) = epoch;
			disp(['Got SNAP data for ',cond,' epoch: ',num2str(epoch.Index)]);
			if ~strcmp(epoch.TrialString,epoch.SNAP.String)
				error('String mismatch'); 
			end
		end
		pipeline.Subjects(subj_i).Clusters(clust_i) = clust;
	end
end

