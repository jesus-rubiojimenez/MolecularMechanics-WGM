function [detrendedData] = wgmDetrend(filename,modeNum,timeSeries,initialTime,finalTime,optPlots)
%% Detrending Whispering Gallery Mode (WGM) sensor data
%
% Jesús Rubio, PostDoc
% University of Exeter
% Created: August 2020
% Last update: June 2021
%
% This function detrends the WGM sensor data, so that we the measurement 
% errors due to factors such as drift and temperature are smoothed out.
%
% For best results, choose (finalTime-initialTime) as a multiple of 20 s,
% and note that the last 10 values are removed to avoid boundary errors.
% 
% Input:
%   - filename (e.g., 'Analyse13Dec19Meas9Ana3.txt')
%   - modeNum: selects sensor tracked mode (1-3)
%       In 7-column files, the first one is time, and the next 6 are
%       3 shifts and 3 FWHM columns.
%   - timeSeries: selects shifts (1) or FWHMs (2)
%   - initialTime: initial time of the relevant part of the trace
%   - finalTime: final time of the relevant part of the trace
%   - optPlots: Detrended plot (2).
%       Choose (1) for plotting the estimated trend and raw data.
%       Choose (0) when no plot is required. 
%
% Output:
%   - detrendedData: structure storing detrended data and the associated times
%
% Remark: to extract maximum information, one should ideally estimate the
% trend and the signal simultaneously. However, it is currently not possible
% to interpret molecular events in a WGM sensor without first detrending the
% traces. In turn, this means that it is not possible to select a suitable
% mathematical representation of the signal unless the raw data is already 
% detrended. 

%% WGM sensor data
rawData=wgmRawData(filename,modeNum,timeSeries,initialTime,finalTime);
timeVar=rawData.rawTime;
depVar=rawData.rawDependentVar;
depVar=depVar*10^6; % From nm to fm

%% Piece-wise smoothing via non-parametric Bayesian estimation with a basis of Gaussian functions
depVar0=depVar(1); % Reference
depVar=depVar-depVar0;

if finalTime-initialTime>20
    windowSize=20; % (in s) Chosen for computational efficiency
else
    windowSize=(finalTime-initialTime);
end

dTime=timeVar(2)-timeVar(1); % Experimental resolution (in s)
movingDetrendNum=floor(windowSize/dTime); % Number of points for the window of local detrending

timeVarLoc=vec2mat(timeVar,movingDetrendNum); % Rows contains the relevant vectors
depVarLoc=vec2mat(depVar,movingDetrendNum);

numPieces=size(timeVarLoc,1); % Number of pieces
locDim=size(timeVarLoc,2); % Dimension of each piece of trace

relFac=10^(-4); % Square of: exp. resolution over width of local fitting Gaussian

i=1:locDim;
j=i;
auxMat=exp(-relFac*(i-j').^2/2);
auxMat=sparse(auxMat);

Q=auxMat*auxMat;
Q=sparse(Q);

invAux=inv(Q+identity(length(Q)));
invAux=sparse(invAux);

tempMat=auxMat+identity(length(auxMat));
tempMat=sparse(tempMat);

finalMat=auxMat*invAux*tempMat;
finalMat=sparse(finalMat);

shiftSmooth=zeros(numPieces,locDim);
detrendedDepVar=zeros(numPieces,locDim);
for k=1:numPieces
    shiftSmooth(k,:)=finalMat*depVarLoc(k,:)';
    detrendedDepVar(k,:)=depVarLoc(k,:)-shiftSmooth(k,:);
end

shiftSmooth=reshape(shiftSmooth',1,movingDetrendNum*numPieces)';
detrendedDepVar=reshape(detrendedDepVar',1,movingDetrendNum*numPieces)';

% We need to remove:
%    - extra zeros origined in how we are splitting the trace
%    - The last 10 values, to avoid boundary errors
timeVar=timeVar(1:length(timeVar)-10,1);
depVar=depVar(1:length(depVar)-10,1);
shiftSmooth=shiftSmooth(1:length(timeVar),1);
detrendedDepVar=detrendedDepVar(1:length(timeVar),1);

%% Detrended data
detrendedData.timeVar=timeVar;
detrendedData.detrendedVar=detrendedDepVar;

%% Plots
if optPlots==0
elseif optPlots==1
    fontsize=30;
    plot(timeVar,depVar,'-','LineWidth',1.5)
    hold on
    plot(timeVar,shiftSmooth,'r','LineWidth',1.5)
    xlim([min(timeVar) max(timeVar)])
    legend('Raw data, $\Delta \lambda^{\mathrm{raw}}$', 'Trend estimate, $\tilde{T}$','Location', 'NorthEast','Interpreter','latex','FontSize',fontsize-10);
    xlabel('$t$ [s]','Interpreter','latex','FontSize',fontsize);
    ylabel('$\Delta \lambda$ [fm]','Interpreter','latex','FontSize',fontsize);
    set(gca,'FontSize',fontsize,'FontName','Times')
    grid on
    hold off
elseif optPlots==2
    fontsize=30;
    plot(timeVar,detrendedDepVar,'LineWidth',1.5)
    xlim([min(timeVar) max(timeVar)])
    xlabel('$t$ [s]','Interpreter','latex','FontSize',fontsize);
    ylabel('$\Delta \lambda$ [fm]','Interpreter','latex','FontSize',fontsize);
    set(gca,'FontSize',fontsize,'FontName','Times')
    grid on
    hold off    
end
end