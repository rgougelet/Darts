clear; close all; clc;
% script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
script_dir = 'G:/darts/';
cd(script_dir);
addpath([script_dir,'deps/'])
out_dir = [script_dir,'mc/'];
mkdir(out_dir);

s = 3;
n_mc_sims = 2000;
n_boot_sims = 2000;
n_obs = 2000;
n_vars = 10;
n_sd = 5;

for sim_i = 1:n_mc_sims
	tic
	svnm = [out_dir,'mc',num2str(sim_i),'.mat'];
	if exist(svnm,'file')
		continue
	end
	runmc(s, n_sd, n_obs, n_vars, n_boot_sims, svnm);
	toc
end

% mxr = mean(arrayfun(@(x) mean(mean((x.XR).^2)), mcs))
% mcs.XR;