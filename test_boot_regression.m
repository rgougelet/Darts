clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
% script_dir = 'G:/darts/';
cd(script_dir);
addpath([script_dir,'deps/'])
% for s = .5:.5:5
% 	out_dir = [script_dir,'mc',num2str(s),'/'];
% 	mkdir(out_dir);
% 	
	n_mc_sims = 10000;
	n_boot_sims = 10000;
	n_obs = 2000;
	n_vars = 10;
	n_sd = 5;
% 	
% 	overwrite = true;
% 	%%
% 	sims_left = 1:n_mc_sims;
% 	if ~overwrite
% 		for sim_i = sims_left
% 			svnm = [out_dir,'mc',num2str(sim_i),'.mat'];
% 			sim_is(sim_i) = ~exist(svnm,'file');
% 		end
% 	else
% 		sim_is = true(size(sims_left));
% 	end
% 	
% 	parfor sim_i = sims_left(sim_is)
% 		tic
% 		svnm = [out_dir,'mc',num2str(sim_i),'.mat'];
% 		runmc(s, n_sd, n_obs, n_vars, n_boot_sims, svnm);
% 		toc
% 	end
% end
%%
% bbhci = nan(n_mc_sims,n_vars+1, 2);
% bhs = nan(n_mc_sims,n_vars+1);
% bhps = nan(n_mc_sims,n_vars+1);
% bbhps = nan(n_mc_sims,n_vars+1);
for s = 1:.5:5
	out_dir = [script_dir,'mc',num2str(s),'/'];
	for sim_i = 1:n_mc_sims
		svnm = [out_dir,'mc',num2str(sim_i),'.mat'];
		if ~exist(svnm,'file');	continue;	end
		mc = load(svnm);
		mcs(sim_i) = mc;
		% 	bhs(sim_i,:) = mc.bh;
		% 	bhps(sim_i,:) = mc.bhp;
		% 	try bc = mc.bhs; catch bc = mc.boot.bhs; end
		% 	% 		bs(sim_i).bbhs = nanmean(bc,1);
		% 	bbhps(sim_i,:) = 2*min([sum(bc>0);sum(bc<0)],[],1)./length(bc);
		% 	% 		bs(sim_i).bbhci = quantile(bc,[.025, .975]);
		% 	bbhci(sim_i,:,:) = quantile(bc,[.025, .975])';
	end
	save([script_dir,'mcs',num2str(s),'.mat'], 'mcs', '-v7.3', '-nocompression')
end

%%
% n_mc_sims = 10000;
% for sim_i = 1:n_mc_sims
% 	corrs(sim_i) = mean(mean(mcs(sim_i).XR(triu(true(size(mcs(sim_i).XR))))));
% % 		r_sqr(sim_i) = mean(mean(mcs(sim_i).)));
% 		r_sqr(sim_i) = mean(mean(mcs(sim_i).XR(triu(true(size(mcs(sim_i).XR))))));
% end
% mean(corrs)
%%
% mean([mcs.adj_r_sqr])
% b = [mcs.b];
% bh = [mcs.bh];
% mean(bh(2:end,:)-b,2)
% mean([-[mcs.bh])

%%
% bhs = cell2mat({bs.bhs})';
% bbhs = cell2mat({bs.bbhs}');
% bhps = cell2mat({bs.bhps})';
% bbhps = cell2mat({bs.bbhps}');
% bbhci = cell2mat({bs.bbhci});

% mxr = mean(arrayfun(@(x) mean(mean((x.XR).^2)), mcs))
% mcs.XR;