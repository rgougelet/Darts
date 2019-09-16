function [EEGout,selnew] = eeg_rejmax(EEG,mx,width)
% reject overload points from continuous data
n = EEG.nbchan;
N = EEG.pnts;

if nargin < 2
    nm = sqrt(sum(EEG.data.^2)/n);
    mx = 3*std(nm);
end
if nargin < 3
    width = EEG.srate;
end

ind = find(max(abs(EEG.data)) > mx);
if isempty(ind)
    EEGout = EEG;
    return
end
sel = zeros(length(ind),2);
for k = 1:length(ind)
    sel(k,:) = [max(1,ind(k)-width),min(EEG.pnts,ind(k)+width)];
end

selnew(1,:) = sel(1,:);
newk = 2;
for k = 2:size(sel,1)
    if sel(k,1) <= selnew(newk-1,2)
        selnew(newk-1,2) = sel(k,2);
    else
        selnew(newk,:) = sel(k,:);
        newk = newk+1;
    end
end    

EEGout = pop_select(EEG,'nopoint',selnew);

