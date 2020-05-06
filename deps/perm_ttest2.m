function [h,p] = perm_ttest2(x,y,varargin)
%% performs exactly as native ttest2 function
% but using permutation statistics, rather
% than parametric...
% bootstrap estimates are returned in the ci
% and stats output, if requested
% this test treats each element within the specified
% dim as independent, e.g. (1,1,:) is not correlated with
% (1,2,:), where dim = 3; if not true then statistical power
% may be reduced.
% can also be used to test medians, as well as means
if nargin > 2
	[varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
	error(message('stats:ttest2:TooFewInputs'));
end

% Process remaining arguments
alpha = 0.05;
tail = 0;    % code for two-sided;
vartype = '';
dim = '';
nperm = '';
stat = '';

if nargin>=3
	if isnumeric(varargin{1})
		% Old syntax
		%    TTEST2(X,Y,ALPHA,TAIL,VARTYPE,DIM)
		alpha = varargin{1};
		if nargin>=4
			tail = varargin{2};
			if nargin>=5
				vartype =  varargin{3};
				if nargin>=6
					dim = varargin{4};
				end
			end
		end
		
	elseif nargin==3
		error(message('stats:ttest2:BadAlpha'));
		
	else
		% Calling sequence with named arguments
		okargs =   {'alpha' 'tail' 'vartype' 'dim' 'nperm' 'stat'};
		defaults = {0.05    'both'    ''      '' 100 'mean'};
		[alpha, tail, vartype, dim, nperm, stat] = ...
			internal.stats.parseArgs(okargs,defaults,varargin{:});
	end
end

if isempty(alpha)
	alpha = 0.05;
elseif ~isscalar(alpha) || alpha <= 0 || alpha >= 1
	error(message('stats:ttest2:BadAlpha'));
end

if isempty(tail)
	tail = 0;
elseif isnumeric(tail) && isscalar(tail) && ismember(tail,[-1 0 1])
	% OK, grandfathered
else
	[~,tail] = internal.stats.getParamVal(tail,{'left','both','right'},'''tail''');
	tail = tail - 2;
end

if isempty(vartype)
	vartype = 1;
elseif isnumeric(vartype) && isscalar(vartype) && ismember(vartype,[1 2])
	% OK, grandfathered
else
	[~,vartype] = internal.stats.getParamVal(vartype,{'equal','unequal'},'''vartype''');
end

if isempty(dim)
	% Figure out which dimension nanmean will work along by looking at x.  y
	% will have be compatible. If x is a scalar, look at y.
	dim = find(size(x) ~= 1, 1);
	if isempty(dim), dim = find(size(y) ~= 1, 1); end
	if isempty(dim), dim = 1; end
	
	% If we haven't been given an explicit dimension, and we have two
	% vectors, then make y the same orientation as x.
	if isvector(x) && isvector(y)
		if dim == 2
			y = y(:)';
		else % dim == 1
			y = y(:);
		end
	end
end

if isempty(nperm)
	nperm = 10000;
end

% Make sure all of x's and y's non-working dimensions are identical.
sizex = size(x); sizex(dim) = 1;
sizey = size(y); sizey(dim) = 1;
if ~isequal(sizex,sizey)
	error(message('stats:ttest2:InputSizeMismatch'));
end

xnans = isnan(x);
if any(xnans(:))
	nx = nansum(~xnans,dim);
else
	nx = size(x,dim); % a scalar, => a scalar call to tinv
end
ynans = isnan(y);
if any(ynans(:))
	ny = nansum(~ynans,dim);
else
	ny = size(y,dim); % a scalar, => a scalar call to tinv
end

s2x = nanvar(x,[],dim);
s2y = nanvar(y,[],dim);
if strcmp(stat,'mean')
	xmean = nanmean(x,dim);
	ymean = nanmean(y,dim);
	difference = xmean - ymean;
elseif strcmp(stat,'median')
	xmedian = nanmedian(x,dim);
	ymedian = nanmedian(y,dim);
	difference = xmedian - ymedian;
end

% Check for rounding issues causing spurious differences
if strcmp(stat,'mean')
	sqrtn = sqrt(nx)+sqrt(ny);
	fix = (difference~=0) & ...                                     % non-zero
		(abs(difference) < sqrtn.*100.*max(eps(xmean),eps(ymean))); % but small
elseif strcmp(stat,'median')
	sqrtn = sqrt(nx)+sqrt(ny);
	fix = (difference~=0) & ...                                     % non-zero
		(abs(difference) < sqrtn.*100.*max(eps(xmedian),eps(ymedian))); % but small
end

if any(fix(:))
	% Fix any columns that are constant, even if computed difference is
	% non-zero but small
	constvalue = min(x,[],dim);
	fix = fix & all(x==constvalue | isnan(x),dim) ...
		& all(y==constvalue | isnan(y),dim);
	difference(fix) = 0;
end

if vartype == 1 % equal variances
	P = cat(dim,x,y);
	n = size(P,dim);
	x_inds = logical([ones(1,nx), zeros(1,ny)]);
	xsubinds = repmat({':'},1,ndims(P));
	ysubinds = repmat({':'},1,ndims(P));
	psubinds = repmat({':'},1,ndims(P));
	disp(['Running permutation t-test using ', num2str(nperm),' iterations']);
	for perm_i = 1:nperm	
		x_inds = x_inds(randperm(n));
		xsubinds{dim} = x_inds;
		ysubinds{dim} = ~x_inds;
		
		px = P(xsubinds{:});
		py = P(ysubinds{:});
		psubinds{dim} = perm_i;
		if strcmp(stat,'mean')
			pxbar(psubinds{:}) = nanmean(px,dim);
			pybar(psubinds{:}) = nanmean(py,dim);
		elseif strcmp(stat,'median')
			pxbar(psubinds{:}) = nanmedian(px,dim);
			pybar(psubinds{:}) = nanmedian(py,dim);
		end
		if any(perm_i == round(nperm*(0:.1:1)))
			disp([num2str(round(100*perm_i/nperm)), '% done'])
		end
	end
	xbar_diffs = pxbar-pybar;
	if tail==0
		p = nansum(abs(xbar_diffs)>=abs(difference),dim)/nperm;
	elseif tail==-1
		p = nansum(xbar_diffs<=difference,dim)/nperm;
	elseif tail==1
		p = nansum(xbar_diffs>=difference,dim)/nperm;
	end
	
elseif vartype == 2 % unequal variances
% 		xbar_diffs = [];
% 		bx = x(:,randi(nx,1,nx*nperm));
% 		by = y(:,randi(ny,1,ny*nperm));
% 		bxbar = nanmean(bx,dim);
% 		bybar = nanmean(by,dim);
% 		xbar_diffs = bxbar-bybar;
% 
% 	if tail==0
% 		p =  nansum(abs(xbar_diffs)>=abs(difference),dim)/nperm;
% 	elseif tail==-1
% 		p = nansum(xbar_diffs<=difference,dim)/nperm;
% 	elseif tail==1
% 		p = nansum(xbar_diffs>=difference,dim)/nperm;
% 	end
% 	
end
h = cast(p <= alpha, 'like', p);
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha

% if vartype == 1 % equal variances
% 	dfe = nx + ny - 2;
% 	se = std(xbar_diffs,1,dim);
% 	ratio = difference ./ se;
% 	if (nargout>3)
% 		stats = struct('tstat', ratio, 'df', cast(dfe,'like',ratio), ...
% 			'sd', se*sqrt(nperm));
% 		if isscalar(dfe) && ~isscalar(ratio)
% 			stats.df = repmat(stats.df,size(ratio));
% 		end
% 	end
% elseif vartype == 2 % unequal variances
% 	dfe = nx + ny - 2;
% 	se_x = std(pxbar,0,dim);
% 	se_y = std(pybar,0,dim);
% 	
% 	if (nargout>3)
% 		stats = struct('tstat', ratio, 'df', cast(dfe,'like',ratio), ...
% 			'sd', cat(dim, s2x, s2y));
% 		if isscalar(dfe) && ~isscalar(ratio)
% 			stats.df = repmat(stats.df,size(ratio));
% 		end
% 	end
% 	if all(se(:) == 0), dfe = 1; end
% end