clear; clc;
scriptdir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis';
cd(scriptdir)

% load XLSX SNAP data
[num,txt,raw] = xlsread('.\behavioral_data.xlsx');
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end
subj_ns = unique(xlsx.subject)';

% trial level accuracy regression
throwtimes = xlsx.throwtime(xlsx.distance < 6 & xlsx.throwtime < 5);
delay = xlsx.delaystilltime(xlsx.distance < 6 & xlsx.throwtime < 5);
pres = xlsx.pres(xlsx.distance < 6 & xlsx.throwtime < 5);
dists = xlsx.distance(xlsx.distance < 6 & xlsx.throwtime < 5);
subjs = categorical(xlsx.subject(xlsx.distance < 6 & xlsx.throwtime < 5));
reg_table = table(throwtimes,delay,dists);
reg_table.pres = categorical(pres);
reg_table.subjs = categorical(subjs);

fit = fitlm(reg_table,'dists~throwtimes+delay+pres+subjs', 'RobustOpts', 'on');

