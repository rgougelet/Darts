function mytopo(data,chanlocs,fignum,dim,compnum)

[nchan,ncomp] = size(data);
if nargin < 4 | dim == 0
    dim = factor2(ncomp);
end
figure(fignum), clf, hold on;

%mn = min(min(data(:,10:end)));
%mx = max(max(data(:,10:end)));

for i = 1:dim(1)
    for j = 1:dim(2)
        k = (i-1)*dim(2)+j;
        if k <= ncomp
            subplot(dim(1),dim(2),k);
            topoplot(data(:,k),chanlocs,'electrodes','off');
%            caxis([mn mx]);
            if nargin >= 5
                title([int2str(k) ' - ' int2str(compnum(k))]);
            else
                title(int2str(k));
            end
        end
    end
end
