function str=uplow(str)
str=lower(str);
idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
str(idx)=upper(str(idx));