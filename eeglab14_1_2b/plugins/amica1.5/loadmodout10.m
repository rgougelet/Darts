function modout = loadmodout10(outdir)
% function mod = loadmodout(outdir)
%
% Load the output of AMICA from outdir.
%
% Example mod structure:
%
% mod = 
%
%    num_models: 3                      # number of models
%      mod_prob: [3x1 double]           # vector of model probabilities 
%             W: [58x58x3 double]       # model weight matrices
%       num_pcs: 58                     # number of principle comps kept
%      data_dim: 58                     # number of channels in data
%     data_mean: [58x1 double]          # data mean
%             S: [58x58 double]         # sphering matrix for all models
%           LLt: [3x1596014 double]     # Likelihood of time pnts by model        
%            LL: [2500x1 double]        # optimization likelihood history
%             c: [58x3 double]          # model centers (means)
%         alpha: [5x58x3 double]        # src density model mix coeff
%            mu: [5x58x3 double]        # src density model means
%         sbeta: [5x58x3 double]        # src density scales
%           rho: [5x58x3 double]        # src density shapes
%            nd: [2500x58x3 double]     # norm of component change in optim
%          svar: [58x3 double]          # model comonent variances
%             A: [58x58x3 double]       # model components (EEG.icawinv)
%       origord: [58x3 double]          # original order of components
%             v: [3x1596014 double]     # posterior model probabilities
%


weights_name = 'W';
sphere_name = 'S';
outdir = [outdir '/'];

fid = fopen([outdir 'gm'],'r');
if fid > 0
    gm = fread(fid,inf,'double');
    num_models = length(gm);
    fclose(fid);
else
    disp('No gm present, setting num_models to 1');
    num_models = 1;
    gm = 1;
end
modout.num_models = num_models;
modout.mod_prob = gm;

fid = fopen([outdir weights_name],'r');
if fid > 0
    W = fread(fid,inf,'double');
    nw2 = length(W)/num_models;
    nw = sqrt(nw2);
    W = reshape(W,nw,nw,num_models);
    fclose(fid);
else
    disp('No W present, exiting');
    return;
end
modout.W = W;
modout.num_pcs = nw;

mnset = 0;
nxset = 0;
fid = fopen([outdir 'mean'],'r');
if fid > 0
    mn = fread(fid,inf,'double');
    nx = length(mn);
    fclose(fid);
    mnset = 1;
    nxset = 1;
else
    disp('No mn present, setting mean to zero');
end


fid = fopen([outdir 'S'],'r');
if fid > 0
    if nxset
        S = fread(fid,[nx nx],'double');
    else
        S = fread(fid,inf,'double');
        nx = sqrt(length(S));
        S = reshape(S,nx,nx);
    end
    fclose(fid);
else
    if nxset
        S = eye(nx);
    else
        disp('No sphere or mn, exiting');
        return;
    end
end
if ~mnset
    mn = zeros(nx,1);
end

modout.data_dim = nx;
modout.data_mean = mn;
modout.S = S;

fid = fopen([outdir 'comp_list'],'r');
if fid > 0
    comp_list = fread(fid,num_models*nw,'int');
    comp_list = reshape(comp_list,nw,num_models);
    fclose(fid);
end
modout.comp_list = comp_list;


fid = fopen([outdir 'LLt'],'r');
if fid > 0
    LLt = fread(fid,inf,'double');
    LLt = reshape(LLt,num_models+1,length(LLt)/(num_models+1));
    fclose(fid);
    LLtset = 1;
else
    disp('LLt not set');
    LLt = 0;
    LLtset = 0;
end
modout.Lht = LLt(1:num_models,:);
modout.Lt = LLt(num_models+1,:);

fid = fopen([outdir 'LL'],'r');
if fid > 0
    LL = fread(fid,inf,'double');
    fclose(fid);
else
    LL = 0;
end
modout.LL = LL;

fid = fopen([outdir 'c'],'r');
if fid > 0
    c = fread(fid,[nw,num_models],'double');
    fclose(fid);
else
    c = zeros(nw,num_models);
end
modout.c = c;

fid = fopen([outdir 'alpha'],'r');
if fid > 0
    alphatmp = fread(fid,inf,'double');
    num_mix = length(alphatmp)/(nw*num_models);
    alphatmp = reshape(alphatmp,num_mix,nw*num_models);
    for h = 1:num_models
        for i = 1:nw
            alpha(:,i,h) = alphatmp(:,comp_list(i,h));
        end
    end

    for h = 1:modout.num_models
        for i = 1:nw
            num_mix_used(i,h) = sum(alpha(:,1,1)>0);
        end
    end
    fclose(fid);
else
    num_mix = 1;
    alpha = ones(num_mix,nw,num_models);
end
modout.alpha = alpha;

fid = fopen([outdir 'mu'],'r');
if fid > 0
    mutmp = fread(fid,num_mix*nw*num_models,'double');
    mutmp = reshape(mutmp,num_mix,nw*num_models);
    for h = 1:num_models
        for i = 1:nw
            mu(:,i,h) = mutmp(:,comp_list(i,h));
        end
    end
    fclose(fid);
else
    mu = zeros(num_mix,nw,num_models);
end
modout.mu = mu;

fid = fopen([outdir 'sbeta'],'r');
if fid > 0
    sbetatmp = fread(fid,num_mix*nw*num_models,'double');
    sbetatmp = reshape(sbetatmp,num_mix,nw*num_models);
    for h = 1:num_models
        for i = 1:nw
            sbeta(:,i,h) = sbetatmp(:,comp_list(i,h));
        end
    end
    fclose(fid);
else
    sbeta = ones(num_mix,nw,num_models);
end
modout.sbeta = sbeta;

fid = fopen([outdir 'rho'],'r');
if fid > 0
    rhotmp = fread(fid,num_mix*nw*num_models,'double');
    rhotmp = reshape(rhotmp,num_mix,nw*num_models);
    for h = 1:num_models
        for i = 1:nw
            rho(:,i,h) = rhotmp(:,comp_list(i,h));
        end
    end
    fclose(fid);
else
    rho = 2*ones(num_mix,nw,num_models);
end
modout.rho = rho;

fid = fopen([outdir 'nd'],'r');
if fid > 0
    nd = fread(fid,inf,'double');
    max_iter = length(nd)/(nw*num_models);
    nd = reshape(nd,max_iter,nw,num_models);
    modout.nd = nd;
    fclose(fid);
    ndset = 1;
else
    modout.nd = 0;
    ndset = 0;
end



num_mod = length(gm);

[modout.mod_prob,gmord] = sort(gm,1,'descend');

W = W(:,:,gmord);
modout.c = c(:,gmord);
alpha = alpha(:,:,gmord);
mu = mu(:,:,gmord);
sbeta = sbeta(:,:,gmord);
rho = rho(:,:,gmord);
modout.Lht = modout.Lht(gmord,:);
modout.comp_list = modout.comp_list(:,gmord);
if ndset
    modout.nd = modout.nd(:,:,gmord);
end

n = modout.num_pcs;

for h = 1:num_mod
    A(:,:,h) = pinv(W(:,:,h)*S(1:nw,:));
end

for h = 1:num_mod
    for i = 1:n
        svar(i,h) = sum( alpha(1:num_mix_used(i,h),i,h) .* (mu(1:num_mix_used(i,h),i,h).^2 + ...
            (gamma(3./rho(1:num_mix_used(i,h),i,h))./gamma(1./rho(1:num_mix_used(i,h),i,h)))./sbeta(1:num_mix_used(i,h),i,h).^2) );
        svar(i,h) = svar(i,h) * norm(A(:,i,h))^2;
    end
    [modout.svar(:,h),origord(:,h)] = sort(svar(:,h),1,'descend');

    modout.A(:,:,h) = A(:,origord(:,h),h);
    modout.W(:,:,h) = W(origord(:,h),:,h);
    modout.alpha(:,:,h) = alpha(:,origord(:,h),h);
    modout.mu(:,:,h) = mu(:,origord(:,h),h);
    modout.sbeta(:,:,h) = sbeta(:,origord(:,h),h);
    modout.rho(:,:,h) = rho(:,origord(:,h),h);
    modout.comp_list(:,h) = modout.comp_list(origord(:,h),h);
    if ndset
        modout.nd(:,:,h) = modout.nd(:,origord(:,h),h);
    end
end
modout.origord = origord;
if LLtset
    for h = 1:num_models
        modout.v(h,:) = 0.4343 * (modout.Lht(h,:) - modout.Lt); % log10 odds of model
    end
end

modout.mod_prob = modout.mod_prob';

for h = 1:num_mod
    for i = 1:nw
        na = norm(modout.A(:,i,h));
        modout.A(:,i,h) = modout.A(:,i,h) / na;
        modout.W(i,:,h) = modout.W(i,:,h) * na;
        modout.mu(:,i,h) = modout.mu(:,i,h) * na;
        modout.sbeta(:,i,h) = modout.sbeta(:,i,h) / na;
    end
end
