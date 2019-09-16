% pop_loadmodout() - import AMICA results
%
% Usage:
%   >> OUTEEG = pop_loadmodout( EEG, filepath, smooth_length, version );
%
% Inputs:
%   EEG            - input dataset
%   filepath       - file path for loadmodout function   
%   smooth_length  - Length of probability smoothing Hanning window (in sec.)
%   version        - AMICA version: "1" for the first version, "2" for
%                    AMICA10
%
% Outputs:
%   OUTEEG         - EEGLAB data structure
%
% Author: Arnaud Delorme & Jason palmer & Ozgur Balkan, SCCN, INC, UCSD

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: pop_loadmodout.m,v $

function [EEG, command] = pop_loadmodout(EEG, filepath, smoothlen, amicaVersion, skipDimensionMatching) 

    if nargin < 1
        help pop_loadmodout;
        return;
    end;
    
    command = '';
    
    if nargin < 2
        % ask user
        filepath = uigetdir('*.*', 'Choose an AMICA output folder -- pop_loadmodout'); 
        if filepath == 0 
            return; 
        end;
    end;
    if nargin<5
        skipDimensionMatching = [];
    end
    
    if 0
        if nargin<4
            editwhat2select = { 'AMICA12' }; %'AMICA first version'};
            uilist = {{'style' 'text' 'string' 'Length of probability smoothing window (sec)'}...
                {'style' 'edit' 'string' '2'}...
                {'style' 'text' 'string' 'AMICA type: '} ...
                {'style' 'popupmenu' 'string' editwhat2select 'value' 1}...
                {'style' 'checkbox' 'string' 'Skip dimension matching of likelihood and current data size, and smoothing' 'value' 0}};
            uigeom = {[1 1] [1 1] 1};
            guititle = 'Load AMICA components -- pop_loadmodout12()';
            result = inputgui(uigeom, uilist, 'pophelp(''pop_loadmodout12'')', guititle, [], 'normal');
            if ~isempty(result)
                smoothlen = eval(['[' result{1} ']']);
                
            else
                if isempty(result)
                    return
                else
                    smoothlen = 2;
                end
            end
            amicaVersion = result{2};
            skipDimensionMatching = result{3};
        else
            
            if isempty(skipDimensionMatching)
                skipDimensionMatching = false;
            end
        end
    end
    
    skipDimensionMatching = true;
    
    %switch amicaVersion
    %    case 1
            mod = loadmodout12([ filepath filesep ]);
            % the variable change in AMICA10
            if isfield(mod,'Lht') && isfield(mod,'Lt') && ~isfield(mod,'LLt')
                mod.LLt = mod.Lht;
                %mod = rmfield(mod,'Lht');
            end
        %case 2
            
        %    mod = loadmodout([ filepath filesep ]);
    %end
    
    EEG.etc.amica = mod;   
    % place AMICA model 1 under regular ICA structure.
    EEG.icaweights = EEG.etc.amica.W(:,:,1);  
    EEG.icasphere  = EEG.etc.amica.S(1:size(EEG.etc.amica.W,1),:); 
    EEG.icawinv    = EEG.etc.amica.A(:,:,1);
    EEG.icaact = [];
    EEG = eeg_checkset(EEG);
    EEG.icaact     = eeg_getica(EEG);
    
    
    if 0 %~skipDimensionMatching
    % recompute probabilities if necessary and smooth probabilities through
    % smoothing LLt.
        if isequal(size(EEG.etc.amica.LLt,2),size(EEG.data,2)) && isequal(size(EEG.etc.amica.LLt,3),size(EEG.data,3))
            fprintf('Smoothing probabilities ... \n')
            [EEG.etc.amica.v_smooth EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlen);
            EEG.etc.amica.v = LLt2v(EEG.etc.amica.LLt);
            
        else
            fprintf('Size of data and probabilities mismatch. Probabilities are being recomputed.. \n')
            EEG.etc.amica.LLt = findLLt(EEG,EEG.etc.amica);
            EEG.etc.amica.v = LLt2v(EEG.etc.amica.LLt);
            fprintf('Smoothing probabilities ... \n')
            [EEG.etc.amica.v_smooth EEG.etc.amica.LLt_smooth] = smooth_amica_prob(EEG.srate,EEG.etc.amica.LLt,smoothlen);
        end
        
        EEG.etc.amica.smooth_length = smoothlen; % store it for further reference.
    end    
    
    %EEG = eeg_checkamica(EEG);
    disp('AMICA components imported successfully, under EEG.etc.amica');
    %command = sprintf('%s = pop_loadmodout(%s,''%s'',%d ,%d, %d)',inputname(1),inputname(1),filepath,smoothlen,amicaVersion, skipDimensionMatching);
    
