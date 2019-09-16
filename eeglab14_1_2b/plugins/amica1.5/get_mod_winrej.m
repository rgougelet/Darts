function wr = get_mod_winrej(LLt,nchan)
% function wr = get_mod_winrej(LLt,nchan)
% Convert LLt into winreg structure with nchan channels for color marking
% of model time points in eegplot
%
% If the active model in LLt is changing too often, this will take a long
% and produce a large winrej structure. It is usually best to "smooth" the
% LLt with a moving average window to eliminate very short duration
% "blips".
%
%    winrej format:    [start end R G B e1 e2 e3 ...] Matrix giving data periods to mark 
%                      for rejection, each row indicating a different period
%                      [start end] = period limits (in frames from beginning of data); 
%                      [R G B] = specifies the marking color; 
%                      [e1 e2 e3 ...] = a (1,nchans) logical [0|1] vector giving 
%                         channels (1) to mark and (0) not mark for rejection.
%
%    Model colors:   1  light_blue
%                    2  light_green
%                    3  light_red
%                    4  light_cyan
%                    5  light_purple
%                    6  light_yellow
%                    7  light_grey
%
%

[M,N] = size(LLt);
[LLr,LLind] = max(LLt);
LLind(LLt(1,:)==0) = 7;

light_blue = [173 216 230] / 255;
light_green = [152 251 152] / 255;
light_red = [255 182 193] / 255;
light_cyan = [209 238 238] / 255;
light_purple = [216 191 216] / 255;
light_yellow = [255 236 139] / 255;
light_grey = [214 214 214] / 255;

colr = {light_blue,light_green,light_red,light_cyan,light_purple,light_yellow,light_grey};

if M > 7
    error('need to add more than 7 colors to the code');
end

ind = 1;
row = 1;
wr = zeros(N,nchan+5);
while ind <= N
    strt = ind;
    c = LLind(ind);    
    
    currcol = colr{c};
    
    % get the endpoint of this segment
    while ind <= N && LLind(ind) == c
        ind = ind + 1;
    end
    stp = ind-1;
    %disp(['ind = ' int2str(ind)]);
    wr(row,1:5) = [strt stp currcol];
    row = row + 1;
end
wr = wr(wr(:,1)>0,:);
