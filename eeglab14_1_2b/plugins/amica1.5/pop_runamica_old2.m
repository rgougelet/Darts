% pop_runamica() - Perform AMICA -- adaptive ICA using multiple models through GUI.

% Optional keywords:
%
%   indir               optional input directory from which to load init
%   num_models          number of models to learn, default = 1
%   num_mix_comps       number of mixture components in source model, def=3
%   share_comps         flag to share components when num_models > 1, def=1
%   comp_thresh         correlation threshold to share component, def=0.97
%   share_start         iteration to start looking for shared components, def=100
%   share_int           number of iterations between sharing checks, def=100
%   numprocs            number or processors (slots) to use, def=8
%   max_iter            maximum number of iterations to perform, def=2000
%   lrate               initial learning rate for natural gradient, def=0.1
%   lratefact           multiplicative factor by which to decrease lrate, def=0.5
%   minlrate            lrate after which to stop, def=1e-8
%   rholrate            initial lrate for shape parameters, def=0.05
%   rho0                initial shape parameter value, def=1.5
%   minrho              minimum shape parameter value, def=1.0
%   maxrho              maximum shape parameter value, def=2.0
%   rholratefact        multiplicative factor by which to dec rholrate, def=0.5
%   do_newton           flag for newton method, default = 1 (do newton)
%   newt_start          for newton method, iter at which to start newton, def=50pop_runamica_old.m
%   newtrate            for newton method, lrate for newton iterations, def=1.0
%   newt_ramp           for newton method, number of iter to ramp up to newtrate, def=10
%   writestep           iteration interval between output writes, def=10
%   write_nd            flag to write history of component update norms, def=1
%   write_llt           flag to write model log likelihoods of time points, def=1
%   do_reject           flag for doing rejection of time points, def=0
%   numrej              for rejection, number of rejections to perform, def=3
%   rejsig              for rejection, number of standard dev of likelihood
%                           below which to reject data
%   rejstart            for rejection, iteration at which to start reject, def=3
%   rejint              for rejection, iteration interval between reject, def=3
%   max_threads         maximum number of threads to use on a node, def=4
%   decwindow           moving average window to detect likelihood decrease, def=1
%   update_A            flag to update mixing matrices, def=1
%   update_c            flag to update model centers, def=1
%   update_gamma        flag to update model probabilities, def=1
%   update_alpha        flag to update source mixture proportions, def=1
%   update_mu           flag to update source mixture mixture locations, def=1
%   update_sbeta        flag to update source mixture scales, def=1
%   invsigmax           maximum value of inverse scale parameters, def=100.0
%   invsigmin           minimum value of inverse scale parameters, def=1e-8
%   do_rho              flag to update shape parameters, def=1
%   load_rej            flag to load LLt to get rejections from, def=0
%   load_param          flag to load parameters, def=0
%   pcakeep             for PCA reduction, number of components to keep, def=chans
%   doscaling           flag to rescale unmixing matrix rows to unit norm, def=1
%   block_size          matrix block size (for block multiplication), def=128
%   qsub                ['on'|'off'] use qsub to run in parallel (default) or run on 
%                       local machine
%
% Disabled:
%   do_mean             flag to remove mean from data, def=1
%   do_sphere           flag to sphere data before ica, def=1
%   doPCA               flag to to PCA dimensionalit reduction, def=0
%   kurt_int            for ext. infomax, iteration interval between calc, def=1
%   kurt_start          for ext. infomax, iter to start kurtosis calc, def=3
%   load_comp_list      flag to load component assignment list, def=0
%   num_kurt            for ext. infomax, number of kurtosis calc, def=5
%   scalestep           iteration interval at which to rescale unmixing rows, def=1
%


function EEG = pop_runamica(EEG)
amicaVersions = {'AMICA12'};
computeVersions = {'In parallel'};
shareCompOptions = {'No' 'Yes'};

callbackAMICAVersion = ['if get(gcbo,''value'') ~= 1,' ...
    'set(findobj(''parent'',gcbf,''tag'',''shareComp''),''enable'',''off'');'...
    'else,' ...
    'set(findobj(''parent'',gcbf,''tag'',''shareComp''),''enable'',''on'');'...' ...
    'end;'];

defaultOutputDirectory = [EEG.filepath filesep 'amicaout'];
uilist = {{'style' 'text' 'string' 'AMICA version'} ...
    {'style' 'popupmenu' 'string' amicaVersions 'value' 1 'callback' callbackAMICAVersion} ...
    {'style' 'text' 'string' 'Output directory'} ...
    {'style' 'edit' 'string' defaultOutputDirectory} ...
    {'style' 'text' 'string' 'Number of Models'} ...
    {'style' 'edit' 'string' '1'} ...
    {'style' 'text' 'string' 'Number of Processors'} ...
    {'style' 'edit' 'string' '4'} ...
    {'style' 'text' 'string' 'Max. number of iterations'} ...
    {'style' 'edit' 'string' '2000'} ...
    {'style' 'text' 'string' 'Share components across models'} ...
    {'style' 'popupmenu' 'string' shareCompOptions 'value' 1 'tag' 'shareComp' 'enable' 'on'} ...
    {'style' 'text' 'string' 'Additional commandline options'} ...
    {'style' 'edit' 'string' ''}};

uigeom = {[1 1] [1 1] [1 1] [1 1] [1 1] [1 1] [1 1]};
guititle = 'Run AMICA -- pop_runamica()';
result = inputgui(uigeom, uilist, 'pophelp(''pop_runamica'')', guititle, [], 'normal');
GUIOptionKeywords = {'outdir' 'num_models' 'numprocs' 'max_iter' 'share_comps'};
if isempty(result)
    return
end

for i = 1:length(GUIOptionKeywords)
    keyword = GUIOptionKeywords{i};
    arglist{2*i-1} = keyword;
    if strcmpi(keyword,'outdir')
        arglist{2*i} = result{2}; 
    end
    if strcmpi(keyword,'num_models')
        arglist{2*i} = str2num(result{3});
    end
    if strcmpi(keyword,'numprocs')
        arglist{2*i} = str2num(result{4});
    end
    if strcmpi(keyword,'max_iter')
        arglist{2*i} = str2num(result{5});
    end
    if strcmpi(keyword,'share_comps')
        arglist{2*i} = fastif(result{6}==1,0,1);
    end
end

additionalOptions = eval( [ '{' result{7} '}' ]);
for i = 1:length(additionalOptions)
    arglist{end+1} = additionalOptions{i};
end

amicaVersion = result{1};
outdir = result{2};

if isfield(EEG,'datfile') && length(EEG.datfile) > 0
    runamica12([EEG.filepath filesep EEG.datfile],'num_chans',EEG.nbchan,arglist{:});
else
    disp('No datfile field found in EEG structure. Writing temp file ...');
    runamica12(EEG.data(:,:),arglist{:});
end        
        
        


if 0 %amicaVersion == 1
    if EEG.trials>1
        %tmpdata = reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3));
        runamica12(EEG.data(:,:),arglist{:});
    else
        runamica12(EEG.data(:,:),arglist{:});
    end
    fprintf('AMICA output is going to be in the folder %s \n',outdir);
else
    % for future AMICA versions
end




