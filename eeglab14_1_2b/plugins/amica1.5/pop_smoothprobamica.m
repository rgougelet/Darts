% changes the smoothed probabilities and likelihoods (v_smooth,LLt_smooth)
function [EEG com] = pop_smoothprobamica(EEG,smoothlength)
com = '';
if ~isfield(EEG.etc,'amica') || ~isfield(EEG.etc.amica,'LLt')
    error('Please check if probabilities are under EEG.etc.amica')
end

if nargin < 2
    prompt = {'Enter the smoothing Hanning window size (in sec): '};
    name = 'Input for probability smoothing for AMICA';
    numlines = 1;
    defaultanswer = {'2'};
    smoothlength = inputdlg(prompt,name,numlines,defaultanswer);
    if isempty(smoothlength) 
        return   
    end 
    
    if isempty(str2num(smoothlength{1}))
        disp('Smoothing window size is not numeric!'); return;
    end
    fprintf('Smoothing probabilities...')
    [EEG.etc.amica.v_smooth,EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,str2num(smoothlength{1}));
    fprintf('Done \n')
    EEG.etc.amica.smooth_length = str2num(smoothlength{1});
    smoothlength = str2num(smoothlength{1});
    
else
    fprintf('Smoothing probabilities...')
    [EEG.etc.amica.v_smooth,EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlength);
    fprintf('Done \n')
    EEG.etc.amica.smooth_length = smoothlength;
end

com = sprintf('%s = pop_smoothamicaprob(%s,%d)',inputname(1),inputname(1),smoothlength);
