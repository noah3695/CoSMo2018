function [] = createCorrGraph(subplotNum, T1, T2, Neuron1, Neuron2, numStim) 
%% Is run to create correlation graphs between Neuron1 and Neuron2
% Run these in the following order to run this function:
% cosmo('CALC_classify'); 
% N1=21; N2=59; %or any 2 neurons you want
% T1=getrow(T,T.neuron==N1);
% T2=getrow(T,T.neuron==N2);
% createCorrGraph(i, T1, T2, N1, N2, numStim)

close all
numTrials = length(T1.spikeNum) / numStim;  

% make vector of indices
numInds = [];
for i=1:numStim
    numInds = [numInds ; i*ones(numTrials,1)];
end

dotSize = 35;
numTrials = length(T1.spikeNum) / numStim;  
%% Use this if you want to plot ALL 5 stimuli
% c = [];
% c(1,:) = [0 171 75] / 255;
% c(2,:) = [0 126 200] / 255;
% c(3,:) = [238 186 13] / 255;
% c(4,:) = [213 82 75] / 255;
% c(5,:) = [109 103 110] / 255;
% colors = [];
% colors = [repmat(c(1,:),numTrials, 1); repmat(c(2,:),numTrials, 1); ...
%     repmat(c(3,:),numTrials, 1); repmat(c(4,:),numTrials, 1); ...
%     repmat(c(5,:),numTrials, 1)];
% scatter(T1.spikeNum,T2.spikeNum, dotSize ,colors,'filled')
% hold on

%% Use this if youw ant to plot only the 2 preferred stimuli
prefDir1 = T1.prefDir(1);
prefDir2 = T2.prefDir(1);
prefDirs = [prefDir1 prefDir2];
numStim = 2;
c = [];
% c(1,:) = [0 171 75] / 255;
% c(2,:) = [0 126 200] / 255;
 c(1,:) = [238 186 13] / 255;
 c(2,:) = [213 82 75] / 255;
% c(5,:) = [109 103 110] / 255;
goodTrials = T1.stimDir == prefDir1;
goodTrials = goodTrials + (T1.stimDir == prefDir2);
T1.thisplot = T1.spikeNum(goodTrials==1);
%c(5,:) = [109 103 110] / 255;
colors = [];
colors = [repmat(c(1,:),numTrials, 1); repmat(c(2,:),numTrials, 1)]; %; ...
    %repmat(c(3,:),numTrials, 1); repmat(c(4,:),numTrials, 1); ...
    %repmat(c(5,:),numTrials, 1)];
    
%% create figure, make it pretty 
scatter(T1.spikeNum(goodTrials==1) ,T2.spikeNum(goodTrials==1), dotSize ,colors,'filled')
hold on


hFig = figure(1);

for k=1:numStim %for each stimulus group
    %# indices of points in this group
    %idx = ( numInds == k );
    idx = ( numInds == prefDirs(k) );
    
    thisSet = [T1.spikeNum(idx == 1), T2.spikeNum(idx == 1)];
    
    %save in Structure
    NeuronCorrAnal{k} = fit_ellipse(thisSet(:,1), thisSet(:,2), hFig, c(k,:)); 
    
    %plot a black cross to show the middle of each population
    p=plot(NeuronCorrAnal{k}.X0_in, NeuronCorrAnal{k}.Y0_in, 'kx', 'MarkerSize', 20,'LineWidth',3);

end

xinfo = xlim;
yinfo = ylim;

%create equal axis for x and y
if xinfo(2) > yinfo(2); yinfo(2) = xinfo(2); else; xinfo(2) = yinfo(2); end;
if xinfo(1) > yinfo(1); yinfo(1) = xinfo(1); else; xinfo(1) = yinfo(1); end;

%create legend before drawPublishAxis artifically
for i = 1:numStim
    LH(i) = plot(nan, nan, 'color', c(i,:));
end
hLeg = legend(LH);
Cond1 = sprintf('Cond %d', prefDir1);
Cond2 = sprintf('Cond %d', prefDir2);
hLeg.String = {Cond1, Cond2};
%hLeg.String = {'Cond 1', 'Cond 2', 'Cond 3', 'Cond 4', 'Cond 5'}; %use this if you want to include all 5 conditions 

[~,hObj]=legend(LH);           % return the handles array
hL=findobj(hObj,'type','line');  % get the lines, not text
N = 4;
set(hL,'linewidth',N)            % set their width property

xLab = sprintf('Neuron %d', Neuron1);
yLab = sprintf('Neuron %d', Neuron2);

drawPublishAxis('xLabel', xLab, ...
    'yLabel', yLab, ...
    'titleStr', 'Correlation Comparison', 'labelFontSize', 20, ...
    'xTick',[xinfo(1) xinfo(2)], 'yTick',[yinfo(1) yinfo(2)]);
keyboard
