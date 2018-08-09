% collect variability

stimDur = 1;
dT = 0.001;
stimRate = 80;
gShared = 0.001;
gIndep = 0;
gAnat = 0.001;
anatOff = -1;

% gShared = 0;
% gIndep = 0;

[gShare] = sharedNoise(gShared,dT,stimDur);
[spkInds,spkVec] = genSpikes(stimDur,stimRate,dT);

numRpts = 5000;
spkCounts = nan(numRpts,1);

for i = 1:numRpts
    if i == 1
        [spksOut] = LIFModel(spkInds,spkVec,gShare,gIndep,gAnat,anatOff,dT,stimDur,1);
    else
        [spksOut] = LIFModel(spkInds,spkVec,gShare,gIndep,gAnat,anatOff,dT,stimDur,0);
    end

spkCounts(i) = numel(spksOut);
end

figure;
hold on;
H = histogram(spkCounts);

counts = H.Values/sum(H.Values);
edges = H.BinEdges(1:end-1)+0.5;
meanCount = mean(spkCounts)
stds = std(spkCounts);

poisson = @(k,lambda) exp(k*log(lambda) - lambda - gammaln(k+1));
normal = @(sigma,k,mu) (1/(2*pi*sigma^2)).*exp(-(k-mu).^2./(2*sigma^2));

countRange = meanCount-3*stds:0.01:meanCount+3*stds;

y_poiss = poisson(countRange,meanCount);
y_normal = normal(stds,countRange,meanCount);

figure;
hold on;
bar(edges,counts);
scalefact = max(counts)/max(y_normal);
scalefact2 = max(counts)/max(y_poiss);
plot(countRange,y_poiss*scalefact2,'b');
plot(countRange,scalefact*y_normal,'r');

