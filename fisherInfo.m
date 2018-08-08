function [FI_corr] = fisherInfo(dataDir,numNeuron,numRpts,numStim)

%% Sort spike counts into mat

% numNeuron = 100;
% numRpts = 50;
% numStim = 5;
plotOn = 1;

countMat = zeros(numNeuron,numRpts,numStim);

for i = 1:numStim
    for j = 1:numNeuron 
        inds = and(neuron==j,stimDir==i);
        countMat(j,:,i) = spikeNum(inds);
    end
end

%% Fit normal distribution to likelihood functions p(spkCount|stim)
% (to get I_fisher for population of independently firing neurons)

normdist = @(x,sigma,mu,scale) (scale)*exp(-((x-mu).^2)./(2*sigma^2));
dnorm_dx = @(x,sigma,mu,scale) ...
    -(scale/(sigma^2)) .* (x-mu) .* exp(-(3.*(x-mu).^2)./(2*sigma^2));

bestPar = nan(numNeuron,3,numStim);

for i = 1:numStim
    for j = 1:numNeuron
        [counts,edges] = histcounts(countMat(j,:,i));
        pCount = counts/sum(counts);
        binSize = (edges(2)-edges(1));
        vals = binSize/2+edges(1):binSize:binSize*numel(pCount)+edges(1);
        
        [bestPar(j,:,i),~] = fminsearch(@fitnormal,[std(vals) mean(vals),max(pCount)],[],vals,pCount);
        
        if i==1 && j==1 && plotOn
            figure;hold on;
            scatter(vals,pCount);
            plot(vals(1):0.01:vals(end),...
                normdist(vals(1):0.01:vals(end),bestPar(j,1,i),bestPar(j,2,i),bestPar(j,3,i)));
            
            plot(vals(1):0.01:vals(end),...
                dnorm_dx(vals(1):0.01:vals(end),bestPar(j,1,i),bestPar(j,2,i),bestPar(j,3,i)));
        end
    end
end

%% Get FI (in style of Zylberberg et al. '16, Hu et al. '14)
% (to get I_fisher for population of correlated neurons)

% load in tuning curves and correlations
load(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)));
load(fullfile(dataDir,sprintf('corr_neuronPairs_%dneurons',numNeuron)));
load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim',numNeuron,numStim)));

stims = 1:numStim;

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
        sig = sigma(j);
        mu = prefDir(j);
        scale = scale;  % for future if scale randomized
        
        dTC_ds(j,i) = dnorm_dx(stims,sig,mu,scale);
    end
end

% Calculate Fisher information for a given stimulus
FI_corr = nan(numStim,1);

for i = 1:numStim
    FI_corr(i) = dTC_ds(:,i)' * corrMat * dTC_ds(:,i);
end


end