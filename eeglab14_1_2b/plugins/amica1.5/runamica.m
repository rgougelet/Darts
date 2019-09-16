% runamica() - Perform AMICA ICA with multiple models and adaptive source
%              densities using asymptotic Newton method.
%              Schedules a run on the cluster using num_procs (keyword, default=8) processors.
% Usage:
%         >> runamica(file,outdir,chans,frames,'Key1',Value1,...);
%
% Inputs:
%
%   file        file name (including path) of file containing floating point data 
%   outdir      name of directory to write output (does not have to exist)
%   chans       number of channels in the data file
%   frames      number of frames in the data file
%
% Optional keywords:
%
%   indir               optional input directory from which to load init
%   num_models          number of models to learn, default = 1
%   num_mix_comps       number of mixture components in source model, def=3
%   pdftype             pdf type for source model, 0=gen. gauss. (default),
%                           1=ext.infomax, 2=gaussian, 3=logistic
%   numprocs            number or processors (slots) to use, def=8
%   max_iter            maximume number of iterations to perform, def=2000
%   lrate               initial learning rate for natural gradient, def=0.1
%   lratefact           multiplicative factor by which to decrease lrate, def=0.5
%   minlrate            lrate after which to stop, def=1e-8
%   rholrate            initial lrate for shape parameters, def=0.05
%   rho0                initial shape parameter value, def=1.5
%   minrho              minimum shape parameter value, def=1.0
%   maxrho              maximum shape parameter value, def=2.0
%   rholratefact        multiplicative factor by which to dec rholrate, def=0.5
%   do_newton           flag for newton method, default = 1 (do newton)
%   newt_start          for newton method, iter at which to start newton, def=50
%   newtrate            for newton method, lrate for newton iterations, def=1.0
%   newt_ramp           for newton method, number of iter to ramp up to newtrate, def=10
%   writestep           iteration interval between output writes, def=10
%   write_nd            flag to write history of component update norms, def=0
%   write_llt           flag to write model log likelihoods of time point, def=1
%   do_reject           flag for doint rejection of time points, def=0
%   numrej              for rejection, number of rejections to perform, def=3
%   rejsig              for rejection, number of standard dev of likelihood
%                           below which to reject data
%   rejstart            for rejection, iteration at which to start reject, def=3
%   rejint              for rejection, iteration interval between reject, def=3
%   kurt_start          for ext. infomax, iter to start kurtosis calc, def=3
%   num_kurt            for ext. infomax, number of kurtosis calc, def=5
%   kurt_int            for ext. infomax, iteration interval between calc, def=1
%   max_threads         maximum number of threads to use on a node, def=4
%   decwindow           moving average window to detect likelihood decrease, def=1
%   update_W            flag to update unmixing matrix, def=1
%   update_c            flag to update model centers, def=1
%   update_gamma        flag to update model probabilities, def=1
%   update_alpha        flag to update source mixture proportions, def=1
%   update_mu           flag to update source mixture mixture locations, def=1
%   update_sbeta        flag to update source mixture scales, def=1
%   invsigmax           maximum value of inverse scale parameters, def=100.0
%   invsigmin           minimum value of inverse scale parameters, def=0.01
%   do_rho              flag to update shape parameters, def=1
%   load_rej            flag to load rejections from, def=0
%   load_param          flag to load parameters, def=0
%   do_mean             flag to remove mean from data, def=1
%   do_sphere           flag to sphere data before ica, def=1
%   doPCA               flag to to PCA dimensionalit reduction, def=0
%   pcakeep             for PCA reduction, number of components to keep, def=chans
%   doscaling           flag to rescale unmixing matrix rows to unit norm, def=1
%   scalestep           iteration interval at which to rescale unmixing rows, def=1
%   block_size          matrix block size (for block multiplication), def=128
%   qsub                ['on'|'off'] use qsub to run in parallel (default) or run on 
%                       local machine
%
% Outputs:
%   
%   To load output use the function loadmodout() after job ends:
%                         
%       mod = loadmodout(outdir);
%
%   mod is a structure containing the output components and density models. mod.A(:,:,h) is the components for model h. 
%   mod.varord(:,h) is the index order of the components in variance order, mod.LLt is the likelihood of time
%   points for each model (if set), mod.LL is the history of the log likelihood over iterations, mod.A(:,:,h)*c(:,h)
%   is the center for model k, mod.W(:,:,h) is the unmixing matrix for model h, and mod.S is the sphering matrix.
%                       
% See also: loadmodout(), getmodLLt(), get_mod_winrej(), plotmodhist()
%
%

function modres = runamica(file,outdir,chans,frames,p1,v1,p2,v2,p3,v3,p4,v4,p5,v5,p6,v6,p7,v7,p8,v8,p9,v9,p10,v10,p11,v11,p12,v12,p13,v13,p14,v14)

if nargin < 4
    help runamica;
    return;
end

modres = [];

%%%%%%%%%%%%%%%%%%%%%%% Set up keyword default values %%%%%%%%%%%%%%%%%%%%%%%%%

epochs = 1;

num_models = 1;
num_mix_comps = 3;
pdftype = 0;
block_size = 128;
max_iter = 2000;
filter_length = 1;
lrate = 0.1;
minlrate = 1e-8;
lratefact = 0.5;
rholrate = 0.05;
rho0 = 1.5;
minrho = 1.0;
maxrho = 2.0;
rholratefact = 0.5;
kurt_start = 3;
num_kurt = 5;
kurt_int = 1;
do_newton = 1;
newt_start = 50;
newtrate = 1.0;
newt_ramp = 10;
do_reject = 0;
numrej = 3;
rejsig = 3.0;
rejstart = 2;
rejint = 3;
max_threads = 12;
writestep = 10;
write_nd = 0;
write_LLt = 1;
decwindow = 1;
maxdecs = 3;
update_W = 1;
update_c = 1;
update_gamma = 1;
update_alpha = 1;
update_mu = 1;
update_sbeta = 1;
invsigmax = 100.0;
invsigmin = 0.00001;
do_rho = 1;
load_rej = 0;
load_W = 0;
load_c = 0;
load_gm = 0;
load_alpha = 0;
load_mu = 0;
load_sbeta = 0;
load_rho = 0;
load_param = 0;
do_mean = 1;
do_sphere = 1;
doPCA = 1;
pcakeep = chans;
pcadb = 30.0;
bytesize = 4;
doscaling = 1;
scalestep = 1;
numprocs = 8;
qsub     = 'on';

%
%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%
%
   for i = 5:2:nargin % for each Keyword
      Keyword = eval(['p',int2str((i-5)/2 +1)]);
      Value = eval(['v',int2str((i-5)/2 +1)]);
      if ~isstr(Keyword)
         fprintf('runamica(): keywords must be strings')
         return
      end
      Keyword = lower(Keyword); % convert upper or mixed case to lower

      if strcmp(Keyword,'num_models')
         if isstr(Value)
            fprintf('runamica(): num_models must be a positive integer');
            return
         else
            num_models = Value;
         end
      elseif strcmp(Keyword,'indir')
         if ~isstr(Value)
            fprintf('runamica(): indir must be a string');
            return
         else
            indir = Value;
            load_param = 1;
            load_W = 1;
            load_c = 1;
            load_gm = 1;
            load_alpha = 1;
            load_mu =1;
            load_sbeta = 1;
            load_rho = 1;
         end
      elseif strcmp(Keyword,'num_mix_comps')
         if isstr(Value)
            fprintf('runamica(): num_mix_comps must be a positive integer');
            return
         else
            num_mix_comps = Value;
         end
      elseif strcmp(Keyword,'pdftype') 
         if isstr(Value)
            fprintf('runamica(): pdftype should be either 0, 1, 2, or 3');
            return
         else
            pdftype = Value;
         end
      elseif strcmp(Keyword,'block_size') 
         if ~isstr(Value)
            fprintf('runamica(): block_size must be a positive integer');
            return
         else 
            block_size = Value;
         end
      elseif strcmp(Keyword,'outdir') 
         if ~isstr(Value)
            fprintf('runamica(): outputdir must be a string');
            return
         else 
            outdir = Value;
         end
      elseif strcmp(Keyword,'max_iter')
         if isstr(Value)
            fprintf('runamica(): max_iter must be a positive integer');
            return
         else
            max_iter = Value;
         end
      elseif strcmp(Keyword,'lrate')
         if isstr(Value)
            fprintf('runamica(): lrate must be a positive number');
            return
         else
            lrate = Value;
         end
      elseif strcmp(Keyword,'minlrate')
         if isstr(Value)
            fprintf('runamica(): minlrate should be lrate to stop at');
            return
         else
            lrate = Value;
         end
      elseif strcmp(Keyword,'lratefact')
         if isstr(Value)
            fprintf('runamica(): lratefact must be a positive number');
            return
         else
            lratefact = Value;
         end
      elseif strcmp(Keyword,'rholrate')
         if isstr(Value)             
            fprintf('runamica(): rholrate must be a number');
            return
         else
            rholrate = Value;
         end
      elseif strcmp(Keyword,'rho0')
         if isstr(Value)
            fprintf('runamica(): rho0 must be a number');
            return
         else
            rho0 = Value;
         end
      elseif strcmp(Keyword,'minrho')
         if isstr(Value)
            fprintf('runamica(): minrho must be a number');
            return
         else
            minrho = Value;
         end
      elseif strcmp(Keyword,'maxrho')
         if isstr(Value)
            fprintf('runamica(): maxrho must be a number');
            return
         else
            maxrho = Value;
         end
      elseif strcmp(Keyword,'rholratefact')
         if isstr(Value)
            fprintf('runamica(): rholratefact must be a number');
         else 
            rholratefact = Value;
         end
      elseif strcmp(Keyword,'invsigmax')
         if isstr(Value)
            fprintf('runamica(): invsigmax must be a number');
         else 
            invsigmax = Value;
         end
      elseif strcmp(Keyword,'invsigmin')
         if isstr(Value)
            fprintf('runamica(): invsigmin must be a number');
         else 
            invsigmin = Value;
         end
      elseif strcmp(Keyword,'do_newton')
         if isstr(Value)
            fprintf('runamica(): do_newton should be 0 or 1');
            return
         else             
            do_newton = Value;
         end
      elseif strcmp(Keyword,'update_w')
         if isstr(Value)
            fprintf('runamica(): update_w should be 0 or 1');
            return
         else             
            update_W = Value;
         end
      elseif strcmp(Keyword,'update_c')
         if isstr(Value)
            fprintf('runamica(): update_c should be 0 or 1');
            return
         else             
            update_c = Value;
         end
      elseif strcmp(Keyword,'update_gamma')
         if isstr(Value)
            fprintf('runamica(): update_gamma should be 0 or 1');
            return
         else             
            update_gamma = Value;
         end
      elseif strcmp(Keyword,'update_alpha')
         if isstr(Value)
            fprintf('runamica(): update_alpha should be 0 or 1');
            return
         else             
            update_alpha = Value;
         end
      elseif strcmp(Keyword,'update_mu')
         if isstr(Value)
            fprintf('runamica(): update_mu should be 0 or 1');
            return
         else             
            update_mu = Value;
         end
      elseif strcmp(Keyword,'update_sbeta')
         if isstr(Value)
            fprintf('runamica(): update_sbeta should be 0 or 1');
            return
         else             
            update_sbeta = Value;
         end
      elseif strcmp(Keyword,'do_rho')
         if isstr(Value)
            fprintf('runamica(): do_rho should be 0 or 1');
            return
         else             
            do_rho = Value;
         end
      elseif strcmp(Keyword,'load_param')
         if isstr(Value)
            fprintf('runamica(): load_param should be 0 or 1');
            return
         else             
            load_W = 1;
            load_c = 1;
            load_gm = 1;
            load_alpha = 1;
            load_mu =1;
            load_sbeta = 1;
            load_rho = 1;        
            load_param = 1;
         end         
      elseif strcmp(Keyword,'write_nd')
         if isstr(Value)
            fprintf('runamica(): write_nd should be 0 or 1');
            return
         else             
            write_nd = Value;
         end
      elseif strcmp(Keyword,'write_llt')
         if isstr(Value)
            fprintf('runamica(): write_llt should be 0 or 1');
            return
         else             
            write_LLt = Value;
         end
      elseif strcmp(Keyword,'newtrate')
         if isstr(Value)
            fprintf('runamica(): newtrate must be a number');
            return
         else
            newtrate = Value;
         end
      elseif strcmp(Keyword,'newt_start')
         if isstr(Value)
            fprintf('runamica(): newt_start must be a number');
            return
         else
            newt_start = Value;
         end
      elseif strcmp(Keyword,'newt_ramp')
         if isstr(Value)
            fprintf('runamica(): newt_ramp must be a number');
            return
         else
            newt_ramp = Value;
         end
      elseif strcmp(Keyword,'do_reject')
         if isstr(Value)
            fprintf('runamica(): do_reject must be 0 or 1');
            return
         else
            do_reject = Value;
         end
      elseif strcmp(Keyword,'numrej')
         if isstr(Value)
            fprintf('runamica(): numrej must be a number');
            return
         else
            numrej = Value;
         end
      elseif strcmp(Keyword,'rejsig')
         if isstr(Value)
            fprintf('runamica(): rejsig must be a number');
            return
         else
            rejsig = Value;
         end
      elseif strcmp(Keyword,'rejstart')
         if isstr(Value)
            fprintf('runamica(): rejstart must be a number');
            return
         else
            rejstart = Value;
         end
      elseif strcmp(Keyword,'rejsig')
         if isstr(Value)
            fprintf('runamica(): rejsig must be a number');
            return
         else
            rejsig = Value;
         end
      elseif strcmp(Keyword,'rejint')
         if isstr(Value)
            fprintf('runamica(): rejint must be a number');
            return
         else
            rejint = Value;
         end
      elseif strcmp(Keyword,'max_threads')
         if isstr(Value)
            fprintf('runamica(): max_threads must be a number');
            return
         else
            maxthreads = Value;
         end
      elseif strcmp(Keyword,'qsub')
         if ~isstr(Value)
            fprintf('runamica(): qsub must be ''on'' or ''off''');
            return
         else
            qsub = Value;
         end
      elseif strcmp(Keyword,'writestep')
         if isstr(Value)
            fprintf('runamica(): writestep must be a number');
            return
         else
            writestep = Value;
         end         
      elseif strcmp(Keyword,'decwindow')
         if isstr(Value)
            fprintf('runamica(): decwindow must be a number');
            return
         else
            decwindow = Value;
         end         
      elseif strcmp(Keyword,'maxdecs')
         if isstr(Value)
            fprintf('runamica(): maxdecs must be a number');
            return
         else
            maxdecs = Value;
         end         
      elseif strcmp(Keyword,'numprocs')
         if isstr(Value)
            fprintf('runamica(): numprocs must be a number');
            return
         else
            numprocs = Value;
         end         
      elseif strcmp(Keyword,'pcakeep')
         if isstr(Value)
            fprintf('runamica(): pcakeep must be a number');
            return
         else
            pcakeep = Value;
         end         
      elseif strcmp(Keyword,'bytesize')
         if isstr(Value)
            fprintf('runamica(): bytesize must be a number');
            return
         else
            bytesize = Value;
         end         
      elseif strcmp(Keyword,'doscaling')
         if isstr(Value)
            fprintf('runamica(): doscaling must be a number');
            return
         else
            doscaling = Value;
         end         
      else
         fprintf(['runamica(): unknown flag: ' Keyword])
         return
      end
   end
%
%%%%%%%%%%%%%%%%%%%%% Check keyword values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%


%%%%%%%%%%%%%%%%%%%%% create the param file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if outdir(end) ~= filesep
    outdir(end+1) = filesep;
end
system(['mkdir ' outdir ' >& /dev/null']);

if isnumeric(file)
    filename = ['tmpdata' num2str(round(rand(1)*100000)) '.fdt' ];
    data2write = file;
    file = fullfile(outdir,filename);
else 
    data2write = [];
end;             

% write data on disk
if ~isempty(data2write)
    disp(sprintf('Writing data file:', filename));
    floatwrite(data2write, fullfile(outdir, filename));
end;             

fid = fopen([outdir 'input.param'],'w+');
if fid < 1
    errordlg('cannot create file in outdir','Bad Input','modal');
    return;
end

fprintf(fid,'files %s\n',file);
fprintf(fid,'outdir %s\n',outdir);
fprintf(fid,'num_models %d\n',num_models);
fprintf(fid,'num_mix_comps %d\n',num_mix_comps);
fprintf(fid,'pdftype %d\n',pdftype);
fprintf(fid,'block_size %d\n',block_size);
fprintf(fid,'max_iter %d\n',max_iter);
fprintf(fid,'num_samples %d\n',epochs);
fprintf(fid,'data_dim %d\n',chans);
fprintf(fid,'field_dim %d\n',frames);
fprintf(fid,'field_blocksize 1\n');
fprintf(fid,'filter_length %d\n',filter_length);
fprintf(fid,'lrate %f\n', lrate);
fprintf(fid,'minlrate %e\n', minlrate);
fprintf(fid,'lratefact %f\n', lratefact);
fprintf(fid,'rholrate %f\n', rholrate);
fprintf(fid,'rho0 %f\n', rho0);
fprintf(fid,'minrho %f\n', minrho);
fprintf(fid,'maxrho %f\n', maxrho);
fprintf(fid,'rholratefact %f\n',rholratefact);
fprintf(fid,'kurt_start %d\n',kurt_start);
fprintf(fid,'num_kurt %d\n',num_kurt);
fprintf(fid,'kurt_int %d\n',kurt_int);
fprintf(fid,'do_newton %d\n',do_newton);
fprintf(fid,'newt_start %d\n',newt_start);
fprintf(fid,'newt_ramp %d\n',newt_ramp);
fprintf(fid,'newtrate %f\n', newtrate);
fprintf(fid,'do_reject %d\n',do_reject);
fprintf(fid,'numrej %d\n',numrej);
fprintf(fid,'rejsig %f\n',rejsig);
fprintf(fid,'rejstart %d\n',rejstart);
fprintf(fid,'rejint %d\n',rejint);
fprintf(fid,'max_threads %d\n',max_threads);
fprintf(fid,'writestep %d\n',writestep);
fprintf(fid,'write_nd %d\n',write_nd);
fprintf(fid,'write_LLt %d\n',write_LLt);
fprintf(fid,'decwindow %d\n',decwindow);
fprintf(fid,'max_decs %d\n',maxdecs);
fprintf(fid,'updateW %d\n',update_W);
fprintf(fid,'update_c %d\n',update_c);
fprintf(fid,'update_gm %d\n',update_gamma);
fprintf(fid,'update_alpha %d\n',update_alpha);
fprintf(fid,'update_mu %d\n',update_mu);
fprintf(fid,'update_beta %d\n',update_sbeta);
fprintf(fid,'invsigmax %f\n',invsigmax);
fprintf(fid,'invsigmin %f\n',invsigmin);
fprintf(fid,'do_rho %d\n',do_rho);
fprintf(fid,'load_rej %d\n',load_rej);
fprintf(fid,'load_W %d\n',load_W);
fprintf(fid,'load_c %d\n',load_c);
fprintf(fid,'load_gm %d\n',load_gm);
fprintf(fid,'load_alpha %d\n',load_alpha);
fprintf(fid,'load_mu %d\n',load_mu);
fprintf(fid,'load_beta %d\n',load_sbeta);
fprintf(fid,'load_rho %d\n',load_rho);
if load_param == 1
    fprintf(fid,'indir %s\n',indir);
end
fprintf(fid,'do_mean %d\n',do_mean);
fprintf(fid,'do_sphere %d\n',do_sphere);
fprintf(fid,'doPCA %d\n',doPCA);
fprintf(fid,'pcakeep %d\n',pcakeep);
fprintf(fid,'pcadb %f\n',pcadb);
fprintf(fid,'byte_size %d\n',bytesize);
fprintf(fid,'doscaling %d\n',doscaling);
fprintf(fid,'scalestep %d\n',scalestep);
fclose(fid);


if strcmpi(qsub, 'on')
    % create the qsub file
    fid = fopen([outdir 'qsub.sh'],'w+');
    
    fprintf(fid,'#!/bin/bash\n');
    fprintf(fid,'#$ -cwd\n');
    fprintf(fid,'#$ -j y\n');
    fprintf(fid,'#$ -S /bin/bash\n');
    fprintf(fid,'#$ -pe mpich %d\n',numprocs);
    
    fprintf(fid,'MPI_DIR=/home/jason/mpich-install\n');
    
    fprintf(fid,...
            ['$MPI_DIR/bin/mpirun -np $NSLOTS -machinefile $TMP/machines /data/common/amica/amica ' ...
             '%s' 'input.param'],outdir);
    fclose(fid);
    
    system(['ssh juggling -n qsub ' outdir 'qsub.sh > ' outdir 'lastqsubid']);
    
    % get the qsub job id
    fid = fopen([outdir 'lastqsubid'],'r');
    str = fscanf(fid,'%s');
    fclose(fid);
    qsubid = sscanf(str(8:end),'%d');
    
    disp(['qsub id = ' int2str(qsubid)]);

else
    system(['/data/common/amica/amica ' fullfile(outdir,'input.param') ]);
    if exist(fullfile(outdir, 'W'))
        disp('Now loading results...');
        disp('WARNING: if you are re-using the same output folder and the function crashed,')
        disp('         it may be past results');
        modres = loadmodout(outdir);
    end;
end;

