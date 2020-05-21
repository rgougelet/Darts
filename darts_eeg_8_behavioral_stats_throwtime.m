
%% todo: verify/remove sqrts
% dv distance = idv delay, idv condition, dv throwtime, i/dv theta, i/dv alpha, (1+subject|throwtime)
% dv throwtime = idv delay, idv condition, i/dv theta, i/dv alpha, (1+subject)
% dv distance = 1 x 2 x 1 x 1 x 1
% dv throwtime = 1 x 2 x 1 x 1

% done, distance by subject per condition
% check theta variability, distance vs theta per condition
% check alpha variability, distance vs alpha per condition





%% STAT - sig - single-trial distance X condition
r = struct2table(load_r.r) ;
ta = r.distance(r.pres == 0);
tp = r.distance(r.pres == 1);
[h,p] = ttest2(ta,tp)
[h,p] = perm_ttest2(ta,tp,'nperm', 10000)

%% PLOT - single-trial distance X condition
clc; close all;
r = struct2table(load_r.r) ;
fig = figure('Position', [0 0 1366 786], 'Visible', 'on'); 
h1 = histogram(r.distance(~~r.pres,:), 'BinEdges',0:1:49, 'FaceAlpha',0.15); hold on;
h2 = histogram(r.distance(~r.pres,:), 'BinEdges',0:1:49, 'FaceAlpha',0.15);
h1v = h1.Values;
h2v = h2.Values; clf;
bar(.5:1:49,[h1v;h2v]', 1,'stacked', 'FaceAlpha',0.15);
ylabel('Trials');
xlabel('Distance (cm)');
xline(mean(r.distance(~~r.pres,:)), 'LineWidth',2, 'Color','b');
xline(mean(r.distance(~r.pres,:)), 'LineWidth',2, 'Color','r');
set(gca,'FontSize',12);
t = title('Dart''s Distance from Bull''s-Eye by Condition');
l = legend({'Target Present','Target Absent', 'Target Present Mean', 'Target Absent Mean'});
set(l,'position',[0.75 0.8 0.1 0.1]);
tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% STAT - sig - single-trial throwtime X condition
r = struct2table(load_r.r) ;
ta = r.throwtime(r.pres == 0);
tp = r.throwtime(r.pres == 1);
[h,p] = ttest2(ta,tp);
[h,p] = perm_ttest2(ta,tp,'nperm', 10000)

%% PLOT - single-trial throwtime X condition
clc; close all;
r = struct2table(load_r.r) ;
fig = figure('Position', [0 0 1366 786], 'Visible','on');
h1 = histogram(r.throwtime(~~r.pres,:), 'BinEdges',0.8:.2:3.6, 'FaceAlpha',0.15); hold on;
h2 = histogram(r.throwtime(~r.pres,:), 'BinEdges',0.8:.2:3.6, 'FaceAlpha',0.15);
h1v = h1.Values;
h2v = h2.Values; clf;
bar(.9:.2:3.6,[h1v;h2v]', 1,'stacked', 'FaceAlpha',0.15);
ylabel('Trials');
xlabel('Throw time (seconds)');
xline(mean(r.throwtime(~~r.pres,:)), 'LineWidth',2, 'Color','b');
xline(mean(r.throwtime(~r.pres,:)), 'LineWidth',2, 'Color','r');
set(gca,'FontSize',12);
t = title('Time Taken to Throw by Condition');
l = legend({'Target Present','Target Absent', 'Target Present Mean', 'Target Absent Mean'});
set(l,'position',[0.75 0.8 0.1 0.1])
xlim([0.75, 3.5]);
tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% STAT - sig - summary statistics distance X condition
r = struct2table(load_r.r) ;
abs_subj_means = []; pres_subj_means = [];
for subj_i = 1:length(subjs_to_include)
	abs_subj_means(end+1,:) = mean(r.distance(r.subject == subj_i & r.pres == 0));
	pres_subj_means(end+1,:) = mean(r.distance(r.subject == subj_i & r.pres == 1));
end
[~,p,~,stats] = ttest(abs_subj_means,pres_subj_means)
[~,p,~,se] = perm_ttest(abs_subj_means,pres_subj_means)
[~,p] = ttest2(abs_subj_means,pres_subj_means)
[~,p] = perm_ttest2(abs_subj_means,pres_subj_means,'nperm', 10000)

%% STAT - sig - summary statistics throwtime X condition
r = struct2table(load_r.r) ;
abs_subj_means = []; pres_subj_means = [];
for subj_i = 1:length(subjs_to_include)
	abs_subj_means(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.pres == 0));
	pres_subj_means(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.pres == 1));
end
[h,p,~,stats] = ttest(abs_subj_means,pres_subj_means)
[xbar,p,~,se] = perm_ttest(abs_subj_means,pres_subj_means)
stats.sd/sqrt(length(abs_subj_means-pres_subj_means))
[h,p] = ttest2(abs_subj_means,pres_subj_means)
[h,p] = perm_ttest2(abs_subj_means,pres_subj_means,'nperm', 10000)

%% distance vs delay correlation
%% plot single-trial distance vs delay correlation by condition
close all;
fig = figure('Position', [0 0 1366 786], 'Visible','off');
t = annotation('textbox',[.5,.98,0,0],...
	'string','Dart''s Distance from Bull''s-Eye by Delay Time',...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);
r = struct2table(load_r.r) ;
abs_delay = r.delay(r.pres==0);
pres_delay = r.delay(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);

% single trial plot
subplot(2,1,1)
plot(abs_delay, sqrt(abs_distance),'b.'); hold on;
title('Target Absent');
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (sqrt(cm))')
l_abs = lsline();
l_abs.Color = 'b';
set(gca,'FontSize',12)

subplot(2,1,2)
plot(pres_delay, sqrt(pres_distance),'r.'); hold on;
title('Target Present');
xlim([2.5, 9.5])
ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (sqrt(cm))')
l_pres = lsline();
l_pres.Color = 'r';
set(gca,'FontSize',12)
tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% near sig - statistics single-trial distance vs delay correlation by condition
r = struct2table(load_r.r) ;
abs_delay = r.delay(r.pres==0);
pres_delay = r.delay(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);
disp('Single Trial')
disp('Target Absent')
[rho, p] = corr(abs_delay,abs_distance); [rho, p]
[rho, p] = corr(abs_delay,abs_distance, 'type','Spearman'); [rho, p]
[rho, p] = corr(abs_delay,abs_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(abs_delay,abs_distance); [rho, p]
[rho, p,~,~] = studentized_perm_corr(abs_delay,abs_distance); [rho, p]
disp('Target Present')
[rho, p] = corr(pres_delay,pres_distance); [rho, p]
[rho, p] = corr(pres_delay,pres_distance, 'type','Spearman'); [rho, p]
[rho, p] = corr(pres_delay,pres_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(pres_delay,pres_distance); [rho, p]
[rho, p,~,~] = studentized_perm_corr(pres_delay,pres_distance); [rho, p]
r1 = corr(abs_delay,abs_distance);
n1 = length(abs_delay);
r2 = corr(pres_delay,pres_distance);
n2 = length(pres_delay);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

%% sig - summary statistic distance vs delay correlation using subject x condition means
r = struct2table(load_r.r) ;
abs_subj_means = []; abs_subj_delays = [];
pres_subj_means = []; pres_subj_delays = [];
for subj_i = 1:length(subjs_to_include)
	for delay = 3:9
		abs_subj_delays(end+1,:) = delay;
		abs_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		rm_abs_subj_means(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		pres_subj_delays(end+1,:) = delay;
		pres_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==1));
		rm_pres_subj_means(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==1));
	end
end
disp('Summary Statistics')
disp('Target Absent')
[rho, p] = corr(abs_subj_delays,abs_subj_means); [rho, p]
[rho, p] = corr(abs_subj_delays,abs_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr(abs_subj_delays,abs_subj_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr(abs_subj_delays,abs_subj_means); [rho, p]
[rho,p] = studentized_perm_corr(abs_subj_delays,abs_subj_means); [rho, p]

disp('Target Present')
[rho, p] = corr(pres_subj_delays,pres_subj_means); [rho, p]
[rho, p] = corr(pres_subj_delays,pres_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr(pres_subj_delays,pres_subj_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr(pres_subj_delays,pres_subj_means); [rho, p]
[rho,p] = studentized_perm_corr(pres_subj_delays,pres_subj_means); [rho, p]

r1 = corr(abs_subj_delays,abs_subj_means);
n1 = length(abs_delay);
r2 = corr(pres_subj_delays,pres_subj_means);
n2 = length(pres_delay);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

[rho, p] = corr(pres_subj_delays,abs_subj_means-pres_subj_means); [rho, p]
[rho, p] = corr(pres_subj_delays,abs_subj_means-pres_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr(pres_subj_delays,abs_subj_means-pres_subj_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr(pres_subj_delays,abs_subj_means-pres_subj_means); [rho, p]
[rho,p] = studentized_perm_corr(pres_subj_delays,abs_subj_means-pres_subj_means); [rho, p]

%% near sig, not really - summary statistic distance vs delay correlation using condition means
r = struct2table(load_r.r) ;
abs_subj_means = []; abs_subj_delays = [];
pres_subj_means = []; pres_subj_delays = [];
for subj_i = 1:length(subjs_to_include)
	for delay = 3:9
		abs_subj_delays(end+1,:) = delay;
		abs_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		rm_abs_subj_means(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==0));
		pres_subj_delays(end+1,:) = delay;
		pres_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==1));
		rm_pres_subj_means(subj_i,delay-2) = mean(r.distance(r.subject==subj_i & r.delay==delay & r.pres==1));
	end
end
abs_cond_subj_means = [];
pres_cond_subj_means = [];
for delay = 3:9
	abs_cond_subj_means(end+1,:) = mean(abs_subj_means(abs_subj_delays==delay));
	pres_cond_subj_means(end+1,:) = mean(pres_subj_means(pres_subj_delays==delay));
end
disp('Summary Statistics')
disp('Target Absent')
[rho, p] = corr((3:9)',abs_cond_subj_means); [rho, p]
[rho, p] = corr((3:9)',abs_cond_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',abs_cond_subj_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',abs_cond_subj_means); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',abs_cond_subj_means); [rho, p]
disp('Target Present')
[rho, p] = corr((3:9)',pres_cond_subj_means); [rho, p]
[rho, p] = corr((3:9)',pres_cond_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',pres_cond_subj_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',pres_cond_subj_means); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',pres_cond_subj_means); [rho, p]
disp('Target Absent - Target Present')
d = abs_cond_subj_means-pres_cond_subj_means;
[rho, p] = perm_corr((3:9)',d); [rho, p]
[rho, p] = corr((3:9)',d); [rho, p]
[rho, p] = corr((3:9)',d, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',d, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',d); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',d); [rho, p]

r1 = corr((3:9)',abs_cond_subj_means);
n1 = length(abs_cond_subj_means);
r2 = corr((3:9)',pres_cond_subj_means);
n2 = length(pres_cond_subj_means);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

%% plot summary statistic distance vs delay correlation using condition means
close all; figure;
subplot(2,1,1)
plot(3:9,abs_cond_subj_means,'ko')
xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

subplot(2,1,2)
plot(3:9,pres_cond_subj_means,'ko')
xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

%% distance vs throwtime correlation

%% STAT - sig - statistics single-trial distance vs throwtime correlation by condition
r = struct2table(load_r.r) ;
abs_throwtime = r.throwtime(r.pres==0);
pres_throwtime = r.throwtime(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);
disp('Single Trial')
disp('Target Absent')
[rho, p] = corr(abs_throwtime,abs_distance); [rho, p]
[rho, p] = corr(abs_throwtime,abs_distance, 'type','Spearman'); [rho, p]
[rho, p] = corr(abs_throwtime,abs_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(abs_throwtime,abs_distance); [rho, p]
[rho, p,~,~] = studentized_perm_corr(abs_throwtime,abs_distance); [rho, p]
disp('Target Present')
[rho, p] = corr(pres_throwtime,pres_distance); [rho, p]
[rho, p] = corr(pres_throwtime,pres_distance, 'type','Spearman'); [rho, p]
[rho, p] = corr(pres_throwtime,pres_distance, 'type','Kendall'); [rho, p]
[rho, p,~,~] = perm_corr(pres_throwtime,pres_distance); [rho, p]
[rho, p,~,~] = studentized_perm_corr(pres_throwtime,pres_distance); [rho, p]
r1 = corr(abs_throwtime,abs_distance);
n1 = length(abs_throwtime);
r2 = corr(pres_throwtime,pres_distance);
n2 = length(pres_throwtime);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

%% PLOT - single-trial distance vs throw time correlation by condition
clc; close all;
fig = figure('Position', [0 0 1366 786], 'Visible','on');
t = annotation('textbox',[.5,.98,0,0],...
	'string','Dart''s Distance from Bull''s-Eye by Throw Time',...
	'FitBoxToText','on',...
	'LineStyle','none',...
	'HorizontalAlignment','center',...
	'FontSize',14);
r = struct2table(load_r.r) ;
abs_throwtime = r.throwtime(r.pres==0);
pres_throwtime = r.throwtime(r.pres==1);
abs_distance = r.distance(r.pres==0);
pres_distance = r.distance(r.pres==1);

% single trial plot
subplot(2,1,1)
plot(abs_throwtime, abs_distance,'b.'); hold on;
title('Target Absent');
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
l_abs = lsline();
l_abs.Color = 'b';
set(gca,'FontSize',12)

subplot(2,1,2)
plot(pres_throwtime, pres_distance,'r.'); hold on;
title('Target Present');
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
l_pres = lsline();
l_pres.Color = 'r';
set(gca,'FontSize',12)
tic; print(fig,t.String,'-dsvg','-r300'); toc;

%% STAT - summary statistic distance vs throwtime correlation using subject means - near sig
r = struct2table(load_r.r) ;
abs_subj_means = []; abs_subj_throwtimes = [];
pres_subj_means = []; pres_subj_throwtimes = [];
for subj_i = 1:length(subjs_to_include)
		abs_subj_throwtimes(end+1,:) = mean(r.throwtime(r.subject==subj_i & r.pres==0));
		abs_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.pres==0));
		pres_subj_throwtimes(end+1,:) = mean(r.throwtime(r.subject==subj_i & r.pres==1));
		pres_subj_means(end+1,:) = mean(r.distance(r.subject==subj_i & r.pres==1));
end
disp('Summary Statistics')
disp('Target Absent')
[rho, p] = corr(abs_subj_throwtimes,abs_subj_means); [rho, p]
[rho, p] = corr(abs_subj_throwtimes,abs_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr(abs_subj_throwtimes,abs_subj_means, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(abs_subj_throwtimes,abs_subj_means); [rho, p]
[rho, p] = studentized_perm_corr(abs_subj_throwtimes,abs_subj_means); [rho, p]

disp('Target Present')
[rho, p] = corr(pres_subj_throwtimes,pres_subj_means); [rho, p]
[rho, p] = corr(pres_subj_throwtimes,pres_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr(pres_subj_throwtimes,pres_subj_means, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(pres_subj_throwtimes,pres_subj_means); [rho, p]
[rho, p] = studentized_perm_corr(pres_subj_throwtimes,pres_subj_means); [rho, p]

r1 = corr(abs_subj_throwtimes,abs_subj_means);
n1 = length(abs_subj_throwtimes);
r2 = corr(pres_subj_throwtimes,pres_subj_means);
n2 = length(pres_subj_throwtimes);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

[rho, p] = corr(abs_subj_throwtimes-pres_subj_throwtimes,abs_subj_means-pres_subj_means); [rho, p]
[rho, p] = corr(abs_subj_throwtimes-pres_subj_throwtimes,abs_subj_means-pres_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr(abs_subj_throwtimes-pres_subj_throwtimes,abs_subj_means-pres_subj_means, 'type','Kendall'); [rho, p]
[rho, p] = perm_corr(abs_subj_throwtimes-pres_subj_throwtimes,abs_subj_means-pres_subj_means); [rho, p]
[rho, p] = studentized_perm_corr(abs_subj_throwtimes-pres_subj_throwtimes,abs_subj_means-pres_subj_means); [rho, p]

%% PLOT - summary statistic distance vs throwtime correlation using subject means
close all; figure;
subplot(3,1,1)
plot(abs_subj_throwtimes,abs_subj_means,'r*')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

subplot(3,1,2)
plot(pres_subj_throwtimes,pres_subj_means,'r*')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

subplot(3,1,3)
plot(abs_subj_throwtimes-pres_subj_throwtimes,abs_subj_means-pres_subj_means,'r*')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

%% STAT  - summary statistic distance vs throwtime correlation using subject x condition means
abs_cond_subj_means = [];
pres_cond_subj_means = [];
for subj_i = 1:10
	for delay = 3:9
		abs_cond_subj_throwtimes(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.delay==delay));
		pres_cond_subj_throwtimes(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.delay==delay));
		abs_cond_subj_means(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.delay==delay));
		pres_cond_subj_means(end+1,:) = mean(r.throwtime(r.subject == subj_i & r.delay==delay));
	end
end
disp('Summary Statistics')
disp('Target Absent')
[rho, p] = corr((3:9)',abs_cond_subj_means); [rho, p]
[rho, p] = corr((3:9)',abs_cond_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',abs_cond_subj_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',abs_cond_subj_means); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',abs_cond_subj_means); [rho, p]
disp('Target Present')
[rho, p] = corr((3:9)',pres_cond_subj_means); [rho, p]
[rho, p] = corr((3:9)',pres_cond_subj_means, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',pres_cond_subj_means, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',pres_cond_subj_means); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',pres_cond_subj_means); [rho, p]
disp('Target Absent - Target Present')
d = abs_cond_subj_means-pres_cond_subj_means;
[rho, p] = perm_corr((3:9)',d); [rho, p]
[rho, p] = corr((3:9)',d); [rho, p]
[rho, p] = corr((3:9)',d, 'type','Spearman'); [rho, p]
[rho, p] = corr((3:9)',d, 'type','Kendall'); [rho, p]
[rho,p] = perm_corr((3:9)',d); [rho, p]
[rho,p] = studentized_perm_corr((3:9)',d); [rho, p]

r1 = corr((3:9)',abs_cond_subj_means);
n1 = length(abs_cond_subj_means);
r2 = corr((3:9)',pres_cond_subj_means);
n2 = length(pres_cond_subj_means);
[p,z] = corr_rtest(r1,r2,n1,n2); p(2)

%% PLOT -  summary statistic distance vs throwtime correlation using subject x condition means
subplot(2,1,1)
plot(3:9,abs_cond_subj_means,'ko')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)

subplot(2,1,2)
plot(3:9,pres_cond_subj_means,'ko')
% xlim([2.5, 9.5])
% ylim([0.75, 7])
xlabel('Throw time (seconds)')
ylabel('Distance (cm)')
lsline;
set(gca,'FontSize',12)
