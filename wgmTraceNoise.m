function [sigma0] = wgmTraceNoise(filename,modeNum,timeSeries,initialTime,finalTime)
%% Random noise in Whispering Gallery Mode (WGM) sensor data
% 
% Jesús Rubio, PostDoc
% University of Exeter
% Created: August 2020
% Last update: Jun 2021
% 
% This function quantifies the measurement noise associated with the 
% detrended WGM sensor data.
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
%   - sigma0: structure with: root of the average moving variance, in fm, 
%             which is how we quantify the overall noise, and root of the 
%             minimum moving variance, also in fm, to quantify the prior 
%             uncertainty of the model amplitudes.

%% Detrended WGM sensor data
[detrendedData] = wgmDetrend(filename,modeNum,timeSeries,initialTime,finalTime,0);
timeVar=detrendedData.timeVar; % Time, in s
depVar=detrendedData.detrendedVar; % Dependent variable, in fm

%% Moving variance
windowSize=0.3; % (in s) Typical duration of relevant molecular events
dTime=timeVar(2)-timeVar(1); % Experimental resolution (in s)
movingVarNum=floor(windowSize/dTime); % Number of points for the window of local variance
movingVar=movvar(depVar,movingVarNum);

%% Quantification of noise
sigma0moleSignal=sqrt(mean(movingVar)); % for estimating the signal
sigma0modelAmplitudes=sqrt(min(movingVar)); % as prior info for the amplitudes of the estimation model

sigma0.sigma0moleSignal=sigma0moleSignal;
sigma0.sigma0modelAmplitudes=sigma0modelAmplitudes;
end

