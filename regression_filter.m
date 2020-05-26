% want to use each epoch's baseline+period to generate its
% regression weights. but we need to pad even the whole
% epoch to avoid filtering edge artifacts. first, we need to
% check to see if the result are the same 
% for filter->regression vs regression-> filter
% both operations are linear, so we'd
% expect order not to matter

tic; load subjects; toc
%%
clearvars -except subjects
srate = 512;

% are the corrected filtered data
% identical to the filtered corrected data?

eog =  [subjects(1).Clusters(1).Data.HEOG;subjects(1).Clusters(1).Data.VEOG];
feog = iirsos.bp(eog,srate,[3 8],[2.75,8.25],.1,0);
r = subjects(1).Clusters(1).Data.Raw;
fr = iirsos.bp(r,srate,[3 8],[2.75,8.25],.1,0);
[cfr, ~] = eog_regression(fr,feog);

cr = eog_regression(r,eog);
fcr = iirsos.bp(cr,srate,[3 8],[2.75,8.25],.1,0);

% [cru,~] = eog_regression(r,eog);
% fcru = iirsos.bp(fcr,srate,[3 8],[2.75,8.25],.1,0);
% 
% [crf,~] = eog_regression(r,feog);
% fcrf = iirsos.bp(crf,srate,[3 8],[2.75,8.25],.1,0);

% plot
starti = randi(length(eog)-1000,1,1);
plt_inds = starti:starti+1000;
close all; figure; rjgplot(cfr(plt_inds));hold on; rjgplot(fcr(plt_inds));

% they are slightly off. likely due to edge artifacts that are mostly
% washed out by the length of the signal
% correcting filtered data must be done with filtered eog
% filtering corrected data must be done with unfiltered eog