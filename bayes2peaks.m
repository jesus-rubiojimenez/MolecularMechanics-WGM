function [estBayesValues,estBayesSignal,estBayesData]=bayes2peaks(filename,modeNum,timeSeries,initialTime,finalTime,valleyLoc,localLength,optPlots)

%% Molecular Mechanics Initiative
%
% Jesus Rubio Jimenez, PostDoc
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
%
% Created: July 2020
% Last update: June 2021
%
% Estimation of PGK two-model signals when measured with a WGM sensor
%
% Input:
%   - filename (e.g., 'Analyse13Dec19Meas9Ana3.txt')
%   - modeNum: selects sensor tracked mode (1-3)
%       In 7-column files, the first one is time, and the next 6 are
%       3 shifts and 3 FWHM columns.
%   - timeSeries: selects shifts (1) or FWHMs (2)
%   - initialTime: initial time of the relevant part of the total trace
%   - finalTime: final time of the relevant part of the total trace
%   - valleyLocation: approximate time where the valley of the 2peak structure
%       can be found.
%   - localLength: length of the portion of trace where the two-model
%       signal is located.
%   - optPlots: plots the estimated signal and detrended data (1).
%       Choose (0) when the plot is not required.
%
% Output:
%   - When optPlots=0: estBayesValues
%   - When optPlots=1: estBayesSignal, estBayesData: structures including
%       the estimates for the 6 parameters of the model, their errors, the
%       estimated signals (2 individual peaks + the associated envelope),
%       their time, the original detrended signal with its time, and the
%       estimates the duration of each event.
%
% Notes:
%
% a) Prior info:
%   - Noise: standard deviation
%   - Two main maxima
%   - Two-peak Gaussian model
%
% b) Experimental data
%
%   - Complete experimental trace (detrended version)
%
% c) units: times in s, amplitudes in fm
tic

%% Information

% Noise (global trace)
sigma0noise=wgmTraceNoise(filename,modeNum,timeSeries,initialTime,finalTime);
sigma0=sigma0noise.sigma0moleSignal;
sigma0amplitudes=sigma0noise.sigma0modelAmplitudes;

% Detrended data (substracts noise due to environmental effects)
detrendedData=wgmDetrend(filename,modeNum,timeSeries,max(valleyLoc-10,initialTime),min(valleyLoc+10,finalTime),0);
timeVar=detrendedData.timeVar;
depVar=detrendedData.detrendedVar;

% Location of the double-peak
dTime=timeVar(2)-timeVar(1);
initialTimeLocal=valleyLoc-localLength/2;
finalTimeLocal=valleyLoc+localLength/2;
iTimeIndex=find(timeVar>initialTimeLocal & timeVar<initialTimeLocal+dTime);
fTimeIndex=find(timeVar>finalTimeLocal & timeVar<finalTimeLocal+dTime);

timeVarLoc=timeVar(iTimeIndex:fTimeIndex,1);
depVarLoc=depVar(iTimeIndex:fTimeIndex,1);

% Information for prior assignment
indexValley=find(timeVarLoc>valleyLoc-dTime/2 & timeVarLoc<valleyLoc+dTime/2);

depVarLeft=depVarLoc(1:indexValley,1);
[~, peaksLocLeft]=findpeaks(depVarLeft);
timePeaksLeft=timeVarLoc(peaksLocLeft);
[sortedAmplLeft, ~]=sort(timePeaksLeft(:),'descend');
T1priorGuess=sortedAmplLeft(1); % First peak
A1priorGuess=depVarLoc(T1priorGuess==timeVarLoc);

depVarRight=depVarLoc(indexValley:length(depVarLoc),1);
[~, peaksLocRightAux]=findpeaks(depVarRight);
depVarRightAux=depVarRight(peaksLocRightAux);
peaksLocRight=zeros(1,length(depVarRightAux));
for i=1:length(depVarRightAux)
    peaksLocRight(i)=find(depVarLoc==depVarRightAux(i));
end
timePeaksRight=timeVarLoc(peaksLocRight);
[sortedAmplRight, ~]=sort(timePeaksRight(:),'ascend');
T2priorGuess=sortedAmplRight(1); % Second peak
A2priorGuess=depVarLoc(T2priorGuess==timeVarLoc);

if A1priorGuess<sigma0 && A2priorGuess<sigma0
    warning('The selected event may lie under the noise level.') % this depends on how the noise level is selected
end

%% Estimation of the signal associated with molecular events

% Numerical parameters
discretenum=10^4;

% Prior probability
ApriorWidth=2*sigma0amplitudes; % noise

A1min=A1priorGuess-ApriorWidth/2;
A1max=A1priorGuess+ApriorWidth/2;
A1=linspace(A1min,A1max,discretenum);

A2min=A2priorGuess-ApriorWidth/2;
A2max=A2priorGuess+ApriorWidth/2;
A2=linspace(A2min,A2max,discretenum);

priorA1=ones(1,length(A1));
priorA1=priorA1/trapz(A1,priorA1);
priorA2=ones(1,length(A2));
priorA2=priorA2/trapz(A2,priorA2);

TpriorWidth=dTime; % experimental resolution

T1min=T1priorGuess-dTime/2;
T1max=T1priorGuess+dTime/2;
T1=linspace(T1min,T1max,discretenum);

T2min=T2priorGuess-dTime/2;
T2max=T2priorGuess+dTime/2;
T2=linspace(T2min,T2max,discretenum);

priorT1=ones(1,length(T1));
priorT2=ones(1,length(T2));
priorT1=priorT1/trapz(T1,priorT1);
priorT2=priorT2/trapz(T2,priorT2);

distPeaks=T2priorGuess-T1priorGuess;

auxProblem=3.5; % from geometry of signal model

W1priorGuess=distPeaks/auxProblem;
W1min=max(W1priorGuess-TpriorWidth/(2*auxProblem),0.001);
W1max=min(W1priorGuess+TpriorWidth/(2*auxProblem),distPeaks);
W1=linspace(W1min,W1max,discretenum);

W2priorGuess=distPeaks/auxProblem;
W2min=max(W2priorGuess-TpriorWidth/(2*auxProblem),0.001);
W2max=min(W2priorGuess+TpriorWidth/(2*auxProblem),distPeaks);
W2=linspace(W2min,W2max,discretenum);

priorW1=1./W1; % scale parameter
priorW1=priorW1/trapz(W1,priorW1);
priorW2=1./W2;
priorW2=priorW2/trapz(W2,priorW2);

% Sampling for Monte Carlo integration
tauMC=10^4; 
samplingA1=zeros(1,tauMC);samplingA2=samplingA1;
samplingT1=samplingA1;samplingT2=samplingA1;
samplingW1=samplingA1;samplingW2=samplingA1;

cumulativeA1=cumsum(priorA1*abs(A1(2)-A1(1)));
cumulativeA2=cumsum(priorA2*abs(A2(2)-A2(1)));
cumulativeW1=cumsum(priorW1*abs(W1(2)-W1(1)));
cumulativeW2=cumsum(priorW2*abs(W2(2)-W2(1)));
cumulativeT1=cumsum(priorT1*abs(T1(2)-T1(1)));
cumulativeT2=cumsum(priorT2*abs(T2(2)-T2(1)));

for runs=1:tauMC
    auxiliarA1=cumulativeA1-rand;
    
    for x=1:length(A1)
        if auxiliarA1(x)>0
            indexA1=x;
            break
        end
    end
    samplingA1(runs)=A1(indexA1);
    
    auxiliarA2=cumulativeA2-rand;
    for x=1:length(A2)
        if auxiliarA2(x)>0
            indexA2=x;
            break
        end
    end
    samplingA2(runs)=A2(indexA2);
    
    
    auxiliarW1=cumulativeW1-rand;
    for x=1:length(W1)
        if auxiliarW1(x)>0
            indexW1=x;
            break
        end
    end
    samplingW1(runs)=W1(indexW1);
    
    auxiliarW2=cumulativeW2-rand;
    for x=1:length(W2)
        if auxiliarW2(x)>0
            indexW2=x;
            break
        end
    end
    samplingW2(runs)=W2(indexW2);
    
    auxiliarT1=cumulativeT1-rand;
    for x=1:length(T1)
        if auxiliarT1(x)>0
            indexT1=x;
            break
        end
    end
    samplingT1(runs)=T1(indexT1);
    
    auxiliarT2=cumulativeT2-rand;
    for x=1:length(T2)
        if auxiliarT2(x)>0
            indexT2=x;
            break
        end
    end
    samplingT2(runs)=T2(indexT2);
end

% Estimates for model parameters
temp=0;
tempA1=0;tempA2=0;tempW1=0;tempW2=0;tempT1=0;tempT2=0;tempCycleTime=0;middlePoint=0;
cycleCut=3*sigma0;
for i=1:tauMC
    temp=temp+model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
    
    % Model parameters
    tempA1=tempA1+samplingA1(i)*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
    tempA2=tempA2+samplingA2(i)*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
    tempW1=tempW1+log(samplingW1(i))*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
    tempW2=tempW2+log(samplingW2(i))*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
    tempT1=tempT1+samplingT1(i)*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
    tempT2=tempT2+samplingT2(i)*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
    
    % Event duration (estimation of a function of the parameters)
    if samplingA1(i)>cycleCut
        aux1=sqrt(log(samplingA1(i)/cycleCut)*2);
    else 
        aux1=0; % to cover events close to chosen noise level
    end
    
    if samplingA2(i)>cycleCut
        aux2=sqrt(log(samplingA2(i)/cycleCut)*2);
    else 
        aux2=0;
    end
    
    tempCycleTime=tempCycleTime+(aux2*samplingW2(i)+aux1*samplingW1(i)+samplingT2(i)-samplingT1(i))*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);

    % Middle point (to find interevent times)
    middlePoint=middlePoint+0.5*(samplingT1(i)+samplingT2(i))*model2gauss(timeVarLoc,depVarLoc,samplingA1(i),samplingA2(i),samplingW1(i),samplingW2(i),samplingT1(i),samplingT2(i),sigma0);
end

estA1=tempA1/temp;
estA2=tempA2/temp;

estT1=tempT1/temp;
estT2=tempT2/temp;

estW1=exp(tempW1/temp);
estW2=exp(tempW2/temp);

estEventDur=tempCycleTime/temp;
estMiddle=middlePoint/temp;

disp('Inference completed')

%% Results (summary of model information)
estimates=[estT1,estT2,estA1,estA2,estW1,estW2,estEventDur,estMiddle];
estBayesValues.estValues=estimates';

%% Plots

% Estimated signal
timeEstSignal=linspace(timeVarLoc(1),timeVarLoc(length(timeVarLoc)),3000);
estSignal=estA1.*exp(-(timeEstSignal-estT1).^2/(2.*estW1.^2))+estA2.*exp(-(timeEstSignal-estT2).^2/(2.*estW2.^2));
estSignal1=estA1.*exp(-(timeEstSignal-estT1).^2/(2.*estW1.^2));
estSignal2=estA2.*exp(-(timeEstSignal-estT2).^2/(2.*estW2.^2));

% Further results (useful for plots)
estBayesSignal.timeEstSignal=timeEstSignal';
estBayesSignal.estSignal=estSignal';
estBayesSignal.estSignal1=estSignal1';
estBayesSignal.estSignal2=estSignal2';
estBayesData.timeVarLoc=timeVarLoc';
estBayesData.depVarLoc=depVarLoc';

if optPlots==1
    fontsize=30;
    figure('Name','Characterising molecular movement');
    plot(timeVarLoc,depVarLoc,'-','LineWidth',2)
    hold on
    plot(timeEstSignal, estSignal1,'r','LineWidth',2)
    plot(timeEstSignal, estSignal2,'r','LineWidth',2)
    plot(timeEstSignal, estSignal,'k-.','LineWidth',2)
    plot(timeVarLoc,sigma0*ones(1,length(timeVarLoc)),'--','Color',[0.5,0.5,0.5],'LineWidth',1.25);
    xlim([min(timeVarLoc) max(timeVarLoc)])
    xlabel('$t$ [s]','Interpreter','latex','FontSize',fontsize);
    ylabel('$\Delta \lambda$ [fm]','Interpreter','latex','FontSize',fontsize);
    set(gca,'FontSize',fontsize,'FontName','Times')
    grid on
    hold off
elseif optPlots==0
end
toc
