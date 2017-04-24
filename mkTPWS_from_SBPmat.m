% mkTPWS.m

% CAUTION: This script can produce multi-GB files if you have large numbers 
% of clicks in one directory

% Script takes output from Simone/Marie's click detector and puts it into 
% a format for use in detEdit.
% Expects standard HARP xwav drive directory structure.


% Output:
% A TTPP.m file containing 5 variables: 
%   MTT: An Nx2 vector of detection start and end times, where N is the
%   number of detections
%   MPP: An Nx1 vector of recieved level (RL) amplitudes.
%   MSP: An NxF vector of detection spectra, where F is dictated by the
%   parameters of the fft used to generate the spectra and any  
%   normalization preferences.
%   MSN: 
%   f = An Fx1 frequency vector associated with MSP

clearvars

% Setup variables:
baseDir = 'E:\uswtrTest'; % directory containing de_detector output
outDir = 'E:\uswtrTest\TPWS'; % directory where you want to save your TPWS files
siteName = 'uswtrTest'; % site name, only used to name the output file
ppThresh = 100; % minimum RL in dBpp. If detections have RL below this
% threshold, they will be excluded from the output file. Useful if you have
% an unmanageable number of detections.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the output file exists, if not, make it
if exist(outDir,'dir')~=7
    fprintf('Creating output directory %s',outDir)
    mkdir(outDir)
end

MTT = [];
MPP = [];
MSP = [];
MSN = [];
f = [];
dirSet = dir(baseDir);
[~,outName] = fileparts(outDir);
for itr0 = 1:length(dirSet)
    if ~strcmp(dirSet(itr0).name,'.')&&~strcmp(dirSet(itr0).name,'..')...
        && ~strcmp(dirSet(itr0).name,outName)
    
        inDir = fullfile(baseDir,dirSet(itr0).name);
        fileSet = what(inDir);
        lfs = length(fileSet.mat);
        for itr2 = 1:lfs
            thisFile = fileSet.mat(itr2);
            
            load(char(fullfile(inDir,thisFile)),'-mat','pos','rawStart',...
                'ppSignal','specClickTf','fs','yFilt')
            
            if ~isempty(pos)
                % Prune detections with low received level if needed.
                keepers = find(ppSignal >= ppThresh);
                ppSignal = ppSignal(keepers);
                pos = pos(keepers,:);
                fileStart = datenum(rawStart(1,:)); 
                
                % ATTN: Calculating detection times.
                % This assumes that times in the position vector "pos" 
                % are relative to the start time of the first raw file. 
                posDnum = (pos(:,1)/(60*60*24)) + fileStart +...
                    datenum([2000,0,0,0,0,0]);
                
                % store to vectors:
                MTT = [MTT; posDnum];
                MPP = [MPP; ppSignal];
                MSP = [MSP; specClickTf(keepers,:)];
                MSN = [MSN; yFilt(keepers,:)];

                if isempty(f)
                    f = 0:(fs/(2*1000))/(size(specClickTf,2)-1):(fs/(2*1000));
                end
                clickTimes = [];
                hdr = [];
                specClickTf = [];
                ppSignal = [];
                yFilt = [];
            end
            fprintf('Done with file %d of %d \n',itr2,lfs)
        end
        
        fprintf('Done with directory %d of %d \n',itr0,length(dirSet))
        
        outFileName = sprintf('%s_%s_TPWS1.mat',siteName,dirSet(itr0).name);
        save(fullfile(outDir,outFileName),'MTT','MPP','MSP','MSN','f',...
            '-v7.3')
        MTT = [];
        MPP = [];
        MSP = [];
        MSN = [];
    end
end