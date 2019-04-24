%% Script timing
clear
start = datetime('now');
debugging = 1;

%% Initialize paths
warning('off','MATLAB:rmpath:DirNotFound')
warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:MKDIR:DirectoryExists')
scriptdir = '/data/mobi/Darts/Analysis/Analysis_May-2018_cont';
cd(scriptdir)
rmpath /data/common/matlab/eeglab/
addpath(genpath('./eeglab14_1_2b'))

if debugging == 1
  inpath = './reg_results_32/';
  cd(inpath)
  mats = dir('*.mat');
  cd(scriptdir)
  outpath = './coupled_osc_32/';
  mkdir(outpath)
  srate = 32;
else
  inpath = './reg_results_512/';
  cd(inpath)
  mats = dir('*.mat');
  cd(scriptdir)
  outpath = './coupled_osc_512/';
  mkdir(outpath)
end

%% Reload regression results to find strongest regressor
for mat_i = 1:length(mats)
  prctdone(mat_i,length(mats));
  load([inpath,mats(mat_i).name]);
  subj_id = mats(mat_i).name(1:3);
  n_comps = length(Bs);
  close all
  for comp_i = 1:n_comps
    Y = 	Ys{comp_i};
    X = Xs{comp_i};
    B = Bs{comp_i};
    S = Ss{comp_i};
    E = Es{comp_i};
    CovB = CovBs{comp_i} ;
    LogL = LogLs{comp_i};
    n_props = size(Y,2);
    
    %% Partial Slopes - Standard Errors Bx
    PartialStdErr = diag(sqrt(CovB));
    PartialStdErr = reshape(PartialStdErr,size(B,1),size(B,2));
    
    %% t-statistics (Beta / PartialStdErr)
    X_with_comp = X;
    appf_ts = B./ PartialStdErr; % ampitude, phase1, phase 2, freq
    if comp_i == 1
      appf_ts = [nan(4,4);appf_ts];
      X_with_comp = [nan(length(X_with_comp),4),X_with_comp];
    elseif comp_i == n_comps
      appf_ts = [appf_ts;nan(4,4)];
      X_with_comp = [X_with_comp,nan(length(X_with_comp),4)];
    else
      appf_ts = [appf_ts(1:(4*(comp_i-1)),:);nan(4,4);appf_ts((4*(comp_i-1)+1):end,:)];
      X_with_comp = [X_with_comp(:,1:(4*(comp_i-1))),nan(length(X_with_comp),4),X_with_comp(:,(4*(comp_i-1)+1):end)];
    end
    
    %% Find max t scores for this comp's amp, phas, freq
    % combines two phase measures into one because they are used together
    apf_ts = [appf_ts(:,1),abs(appf_ts(:,2))+abs(appf_ts(:,3)),appf_ts(:,4)];
    [~,max_apf_ts] = maxk(abs(apf_ts),1);
    [~,max_appf_ts] = maxk(abs(appf_ts),1);
    
    %% Label coupling
    to_Hzs(comp_i,:) = comp_i;
    to_props(comp_i,:) = ['A', 'P', 'F']; % fixed values
    from_Hzs(comp_i,:) = ceil(max_apf_ts/n_props);
    from_apf = mod(max_apf_ts,n_props)+1;
    from_props_temp = ['F','A', 'P', 'P'];
    from_props(comp_i,:)  = from_props_temp(from_apf);
    
    %       x = X_with_comp(:,coupled(prop_i));
    %       y = Y(:,prop_i);
    %       plot(x);
    %       hold on;
    %       plot(y);
    %       figure;
    %       scatter(x,y);
    %       [~,p] = corr(x,y);
    %       [~,p] = corr(diff(x),diff(y));
    %       [~,p] = corr(x,y, 'type', 'spearman');
    %       [~,p] = corr(diff(x),diff(y), 'type','spearman');
    %       disp([to_prop_str, from_prop_str, num2str(comp_i), ':', num2str(from_Hz(prop_i))])
    
    %% p-values
    %     pVals = 2*(1-tcdf(abs(tRatios),length(X_with_comp)-2));
    %     pVals(pVals>0.05/numel(pVals))=1;
    %     pVals(pVals<0.05/numel(pVals))=0;
    parsave([outpath, subj_id,'_',num2str(comp_i),'.mat'],...
      {to_Hzs,from_Hzs, to_props,from_props},...
      {'to_Hzs','from_Hzs', 'to_props','from_props'});
  end
end

%% Script timing
stop = datetime('now');
disp(['Start ',datestr(start)])
disp(['End ',datestr(stop)])
disp(['Runtime ', char(duration(stop-start))])