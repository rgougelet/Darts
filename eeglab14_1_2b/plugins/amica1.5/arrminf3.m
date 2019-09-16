function [mo,ord] = arrminf(mi0,maxpass,fignum,ord0)
global n mi ord
mi = mi0;
[m,n] = size(mi);
if nargin < 4
    ord = 1:n;
else
    ord = ord0;
end
minimprov = 1;
still_changing = 1;
pass = 0;
if nargin < 2
    maxpass = 3;
end
if nargin < 3
    f = figure;
else
    f = fignum;
end
figure(f), imagesc(mi);%, colormap hot;


if nargin < 4
sorted = 0;
while (1 & sorted == 0)
    sorted = 1;
    for i = 1:(n-1)
       if sum(mi(:,i)) < sum(mi(:,i+1))
           doswap(i,i+1);
           %tmp = ord(i);
           %ord(i) = ord(i+1);
           %ord(i+1) = tmp;
           sorted = 0;
       end     
    end
    %numpasses = numpasses + 1;
    %disp(['finished pass ' int2str(numpasses)]);
end
end
p0=0.9;
while still_changing & pass < maxpass
    %pass = pass+1;
    %disp(['pass ' int2str(pass)]);
    still_changing = 1;

    rnd = 1:n;%randperm(n);
    for k = rnd
        
        % move in the direction or info c.o.m.
        if k < n && k > 1 && sum(mi(1:k-1,k).*(((k-1):-1:1)').^p0) < sum(mi(k+1:end,k).*((1:n-k)').^p0)            
            doswap(k,k+1);
        elseif k > 1 && k < n && sum(mi(1:k-1,k).*(((k-1):-1:1)').^p0) > sum(mi(k+1:end,k).*((1:n-k)').^p0)
            doswap(k,k-1);
        end
        figure(f), imagesc(mi), colorbar;
        pause(0.3);
    end
end

mo = mi;



function doswap(i,j)
global mi ord
swaprow(i,j);
swapcol(i,j);
tmp = ord(i);
ord(i) = ord(j);
ord(j) = tmp;

function swaprow(i,j)
global mi
tmp = mi(i,:);
mi(i,:) = mi(j,:);
mi(j,:) = tmp;

function swapcol(i,j)
global mi
tmp = mi(:,i);
mi(:,i) = mi(:,j);
mi(:,j) = tmp;