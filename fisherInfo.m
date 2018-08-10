function [FI_corr,dprime] = fisherInfo(dataDir,numNeuron,numRpts,numStim,stims,popType)
% Calculates Fisher information for a population of independently firing
% neurons and for correlated populations
%
% Usage: [FI_corr,dprime] = fisherInfo(dataDir,numNeuron,numRpts,numStim,stims)
%        [FI_corr,dprime] = fisherInfo(dataDir,1000,80,17,[1,1.25,...])


%% Load in tuning curves and correlations
load(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)));
load(fullfile(dataDir,sprintf('corr_neuronPairs_%dneurons_%sPopulation',numNeuron,popType)));
load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron,numStim,popType)));

%% Sort spike counts into mat

plotOn = 1;

countMat = zeros(numNeuron,numRpts,numStim);

for i = 1:numStim
    for j = 1:numNeuron 
        inds = and(neuron==j,stimDir==stims(i));
        countMat(j,:,i) = spikeNum(inds);
    end
end

%% Fit normal distribution to likelihood functions p(spkCount|stim)
% (to get I_fisher for population of independently firing neurons)

normdist = @(x,sigma,mu,scale) (scale)*exp(-((x-mu).^2)./(sigma));
dnorm_dx = @(x,sigma,mu,scale) ...
    -(scale/sigma) .* 2*(x-mu) .* exp(-((x-mu).^2)./sigma);

bestPar = nan(numNeuron,3,numStim);

for i = 1:numStim
    for j = 1:numNeuron
        [counts,edges] = histcounts(countMat(j,:,i),10);
        pCount = counts/sum(counts);
        binSize = (edges(2)-edges(1));
        vals = binSize/2+edges(1):binSize:binSize*numel(pCount)+edges(1);
        
        [bestPar(j,:,i),~] = fminsearch(@fitnormal,[std(vals) mean(vals),max(pCount)],[],vals,pCount);
        
        if i==1 && j==1 && plotOn
            figure;hold on;
            set(gca,'YLim',[0,1],'XLim',[60 120]);
            bar(vals,pCount);
            plot(vals(1):0.01:vals(end),...
                normdist(vals(1):0.01:vals(end),bestPar(j,1,i),bestPar(j,2,i),bestPar(j,3,i)));
            
            plot(vals(1):0.01:vals(end),...
                dnorm_dx(vals(1):0.01:vals(end),bestPar(j,1,i),bestPar(j,2,i),bestPar(j,3,i)));
            text(vals(2),0.7,['sigma:',num2str(bestPar(j,1,i)),...
                ', mu:',num2str(bestPar(j,2,i)),', scale:',num2str(bestPar(j,3,i))]);
        end
    end
end





% keyboard;

%% Get FI (in style of Zylberberg et al. '16, Hu et al. '14)
% (to get I_fisher for population of correlated neurons)

% Sort interneuronal correlations into square mat
corrMat = nan(numNeuron,numNeuron);

for i = 1:numNeuron
    ends = find(neuron2 == numNeuron);
    starts = [1;ends+1];
    corrMat(i,i:end) = corrVar(starts(i):ends(i));
    corrMat(i:end,i) = corrVar(starts(i):ends(i));
end

% Derivative of given neuron's TC about presented stimulus
dTC_ds = nan(numNeuron,numStim);

for i = 1:numStim
    for j = 1:numNeuron
        dTC_ds(j,i) = dnorm_dx(stims(i),sigma(j),prefDir(j),scale(j));
    end
end

% Calculate Fisher information for a given stimulus
FI_corr = nan(numStim,1);

for i = 1:numStim
    FI_corr(i) = dTC_ds(:,i)' * corrMat * dTC_ds(:,i);
end

%% Get dprime for the population for a given pair of stimuli
dStim = stims(2)-stims(1);

dprime = dStim.*sqrt(FI_corr);

figure;
bar(stims,dprime);
ylabel(['d',char(39),' (AU)']);
xlabel('Stimulus (AU)');

%% get entropy of the stimulus

pStim = 1/numStim;

entStim = -sum(pStim*log(pStim));


end