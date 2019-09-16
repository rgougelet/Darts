function ld = getloads(hostlist)
%hostlist = {'computing-0-0','computing-0-1','computing-0-2', 'computing-0-3','computing-0-4','computing-0-5','computing-0-6','computing-0-7','computing-0-8','computing-0-9','computing-0-10'};nhost = length(hostlist);
nhost = length(hostlist);
clear res
for k = 1:nhost
    try
        [a,res{k}] = system(['ssh ' hostlist{k} ' -n cat /proc/loadavg']);
        res2 = res{k}(155:end);
        a = sscanf(res2,'%f');
        ld(k) = a(3);
    catch
        ld(k) = inf;
    end
end