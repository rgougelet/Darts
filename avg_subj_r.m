function r = avg_subj_r(r, var1, var2)
r = r(r.(var1),:);
r = r(r.(var2),:);
vn = r.Properties.VariableNames(2:end);
[G,~] = findgroups(r.Subject);
rarr = splitapply(@mean,table2array(r(:,2:end)),G);
r = array2table(rarr, 'VariableNames',vn);

% r = load_r;
% r = r(r.Cluster==1,:);
% vn = r.Properties.VariableNames(3:end);
% [G,results] = findgroups(table(r.Subject,r.TargetAbsent));
% rarr = table2array(r(:,3:end));
% rres = splitapply(@mean,rarr,G);
% rest = array2table(rres, 'VariableNames',vn);
% [results, rest]
