function [rawData] = wgmRawData(filename,modeNum,timeSeries,initialTime,finalTime)
%% Whispering Gallery Mode (WGM) sensor data
% 
% Jesús Rubio, PostDoc
% University of Exeter
% Created: August 2020
% Last update: August 2020
% 
% This function retrives the WGM sensor data from the .txt files. 


%% How to use this function
%
% [seriesPar,groupEventData,eventData,dataDetrended,rawData]=wgmData(filename,modeNum,timeSeries,sTime,eTime,negFlag,analysisType,durThresh)
% 
% Input:
%   - filename (e.g., 'Analyse13Dec19Meas9Ana3.txt')
%   - modeNum: selects sensor tracked mode (1-3)
%       In 7-column files, the first one is time, and the next 6 are
%       3 shifts and 3 FWHM columns.
%   - timeSeries: selects shifts (1) or FWHMs (2)
%   - initialTime: initial time of the relevant part of the trace 
%   - finalTime: final time of the relevant part of the trace
%
% Output:
%   - rawData: structure storing unprocessed data and the associated times
%  
% Notes: 
%   - Times are measured in s
%   - Raw shifts and FWHMs measured in nm
% 
%% Available data
% 
% Narima's experiment (DNA)
%
% Different concentrations of NaCl
%
%   Analyse17Sep19Meas8Ana0.txt 20ml NaCl
%   Analyse17Sep19Meas9Ana0.txt 50ml NaCl
%   Analyse17Sep19Meas10Ana0.txt 100ml NaCl
%   Analyse17Sep19Meas10Ana1.txt
%   Analyse17Sep19Meas10Ana2.txt
%   Analyse17Sep19Meas11Ana0.txt 150ml NaCl
%   Analyse17Sep19Meas12Ana0.txt 200ml NaCl
%   Analyse17Sep19Meas12Ana1.txt
%
% Different concentrations of DNA molecules
%
%   Analyse17Sep19Meas2Ana0.txt 2nM
%   Analyse17Sep19Meas3Ana0.txt 10nM
%   Analyse17Sep19Meas4Ana0.txt 200nM
%   Analyse17Sep19Meas5Ana0.txt 400nM
%   Analyse17Sep19Meas7Ana0.txt 1um = 1000mM
%
% Initial sample
%
%   Analyse17Sep19Meas1Ana0 
%        Docking DNA binding to gold; measurement starts after 55s
%   Analyse17Sep19Meas6Ana0 
%        Complementary DNA transient interaction 
%
%
% Simona's experiment (PGK)
%
% Double-peak structure
%
%   Analyse13Dec19Meas10AnaX 
%       with X= 0, ..., 21
%       Truly useful files: 3, 9, 21
%   Analyse13Dec19Meas9Ana3
%
% Negative control
% 
%   Analyse06Mar20Meas10Ana0
%   Analyse06Mar20Meas10Ana1
%   Analyse06Mar20Meas10Ana2
%   Analyse06Mar20Meas11Ana0
%   Analyse06Mar20Meas11Ana1

%% Number line before the variable name header
rearrangeInfo=regexp(fileread(filename),'\n','split');
headersRefNumber=find(contains(rearrangeInfo,'End'));

%% Variable name header
varNameHeader=cell2mat(rearrangeInfo(headersRefNumber+1));

%% Complete set of raw data
fid=fopen(filename); 
dataCell = textscan(fid, '%f%f%f', 'CommentStyle', {'File', varNameHeader}, 'CollectOutput',1);
fclose(fid);
rawDataComplete=dataCell{1};

% To handle files with data from more than one EM mode:
if length(rearrangeInfo) < length(rawDataComplete(:,1))
    fid=fopen(filename,'rt');
    dataCell = textscan(fid, '%f%f%f%f%f%f%f', 'CommentStyle', {'File', varNameHeader}, 'CollectOutput',1);
    rawDataComplete=dataCell{1};
    fclose(fid);
end

%% Selected raw data
timeStep=rawDataComplete(2,1)-rawDataComplete(1,1);
iTimeIndex=find(rawDataComplete(:,1)>initialTime & rawDataComplete(:,1)<initialTime+timeStep);
fTimeIndex=find(rawDataComplete(:,1)>finalTime & rawDataComplete(:,1)<finalTime+timeStep);
numColumns=size(rawDataComplete,2);
rawDataSelected=rawDataComplete(iTimeIndex:fTimeIndex,1:numColumns);

%% Final raw data into a structure
rawTime=rawDataSelected(:,1);
rawData.rawTime=rawTime;
if modeNum==1 && numColumns==3
    if timeSeries==1
        rawShift=rawDataSelected(:,2);
        rawData.rawDependentVar=rawShift;
    elseif timeSeries==2
        rawFWHM=rawDataSelected(:,3);
        rawData.rawDependentVar=rawFWHM;
    end
elseif modeNum==1 && numColumns==7
    if timeSeries==1
        rawShift=rawDataSelected(:,2);
        rawData.rawDependentVar=rawShift;
    elseif timeSeries==2
        rawFWHM=rawDataSelected(:,5);
        rawData.rawDependentVar=rawFWHM;
    end
elseif modeNum==2 && numColumns==7
    if timeSeries==1
        rawShift=rawDataSelected(:,3);
        rawData.rawDependentVar=rawShift;
    elseif timeSeries==2
        rawFWHM=rawDataSelected(:,6);
        rawData.rawDependentVar=rawFWHM;
    end
elseif modeNum==3 && numColumns==7
    if timeSeries==1
        rawShift=rawDataSelected(:,4);
        rawData.rawDependentVar=rawShift;
    elseif timeSeries==2
        rawFWHM=rawDataSelected(:,7);
        rawData.rawDependentVar=rawFWHM;
    end
else
    
end