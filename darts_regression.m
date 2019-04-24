%% Script timing
clear
clc
start = datetime('now');
debugging = 1;

%% Initialize paths
warning('off','MATLAB:rmpath:DirNotFound')
warning('off','MATLAB:dispatcher:nameConflict')
scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
cd(scriptdir)
rmpath /data/common/matlab/eeglab/
addpath(genpath('./eeglab14_1_2b'))

if debugging == 1
	inpath = './reg_mats_32/';
	cd(inpath)
	mats = dir('*.mat');
	cd(scriptdir)
	outpath = './reg_results_32/';
	mkdir(outpath)
	srate = 32;
else
	inpath = './reg_mats_512/';
	cd(inpath)
	mats = dir('*.mat');
	cd(scriptdir)
	outpath = './reg_results_512/';
	mkdir(outpath)
end

%% Multivariate multiple regression
for matfile = 1:length(mats)
	
	load([inpath,mats(matfile).name]);
	
	% Reshape X carefully
	[n_props, n_points, n_trials, n_comps] = size(reg_mat);

	x1 = [];
	for trial = 1:n_trials
		x1 = cat(2, x1, squeeze(reg_mat(:,:,trial,:)));
	end
	x2 = [];
	for comp_i = 1:n_comps
		x2 = cat(1, x2, squeeze(x1(:,:,comp_i)));
	end
	
	nans = [];
	for x_i = 1:length(x2)
		nans(end+1) = sum(isnan(x2(:,x_i)));
	end

	% regression model for each component:frequency
	Xs = {};
	Ys = {};
	Bs = {};
	Ss = {};
	Es = {};
	CovBs = {};
	LogLs = {};
	parfor comp_i = 1:n_comps
		prctdone(comp_i-1,n_comps);
		X = x2';
		comp_start_pos= n_props*(comp_i-1)+1;
		comp_end_pos= n_props*(comp_i-1)+n_props;
		Y = X(:,comp_start_pos:comp_end_pos);
		[n,d] = size(Y);
		X(:,comp_start_pos:comp_end_pos) = [];
		
		Ys{comp_i} = Y;
		Xs{comp_i} = X;
		[Bs{comp_i},...
			Ss{comp_i},...
			Es{comp_i} ,...
			CovBs{comp_i} ,...
			LogLs{comp_i}] = mvregress(X,Y);

	% 	% use for planning vs. memory
	% 	chisq = 2*(loglik2-loglik)
	%   p = 1-chi2cdf(chisq, nregions-1)
		
	% 	% plot
	% 	if comp_ind == 1
	% 		B = [zeros(4,4);B];
	% 	elseif comp_ind == n_comps
	% 		B = [B;zeros(4,4)];
	% 	else
	% 		B = [B(1:(4*(comp_ind-1)),:);zeros(4,4);B((4*(comp_ind-1)+1):end,:)];
	% 	end
	% 	Bs = cat(1,Bs,B);
	end
	parsave([outpath, mats(matfile).name],{Ys,Xs,Bs,Ss,Es,CovBs, LogLs},{'Ys','Xs','Bs','Ss','Es','CovBs','LogLs'});
end

%% Correlation
% [R, P] = corrcoef(x);

%% VAR Granger Causality Approach
% [A,SIG] = tsdata_to_var(x,2);
% [G,info] = var_to_autocov(A,SIG);
% time-domain pairwise-conditional causalities
% F = autocov_to_pwcgc(G);

% pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % note arguments!
% sig  = significance(pval,alpha,mhtc);
% 
% % Plot time-domain causal graph, p-values and significance.
% 
% figure(2); clf;
% subplot(1,3,1);
% plot_pw(F);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha)])
% 
% cd = mean(F(~isnan(F)));
% fprintf('\ncausal density = %f\n',cd);

%% Script timing
stop = datetime('now');
disp(['Start ',datestr(start)])
disp(['End ',datestr(stop)])
disp(['Runtime ', char(duration(stop-start))])