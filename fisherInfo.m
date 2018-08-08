function [] = fisherInfo()

% Sort spike counts into mat

numNeuron = 100;
numRpts = 50;
numStim = 5;
plotOn = 1;

countMat = zeros(numNeuron,numRpts,numStim);

for i = 1:numStim
    for j = 1:numNeuron 
        inds = and(neuron==j,stimDir==i);
        countMat(j,:,i) = spikeNum(inds);
    end
end

% Fit normal distribution to likelihood functions

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

% Get FI for each stimulus

stim = 1;

% Measure fisher information
FI = @() 




end