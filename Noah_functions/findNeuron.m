function [T1,T2,N1,N2,numStim] = findNeuron(t, estimate, estimate2)

[f, f2] = sort(abs( t.corrMean-estimate )); %corrMean estimate
corrM_est = f2(1:5);
[f11, f21] = sort(abs( t.corrVar-estimate2 )); %corrMean estimate
corrV_est = f21(1:5);

%find the similar value
notgood = 1;
while notgood == 1
    for j = 1:length(corrM_est)
        for k = 1:length(corrV_est)
            if corrM_est(j) == corrV_est(k)
                notgood = 0;
                N1 = t.neuron1(corrM_est(j));
                N2 = t.neuron2(corrM_est(j));
                                
                baseDir = fileparts(which('cosmo.m'));
                dataDir = fullfile(baseDir,'data');
                numNeuron=100;
                numStim = 5;
                popType='mixture';
                
                T = load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron,numStim,popType)));
                
                T1=getrow(T,T.neuron==N1);
                T2=getrow(T,T.neuron==N2);
                return
            end
        end
    end
    if j == length(corrM_est) && k == length(corrV_est)
        error %stop it because you've done fucked up (haven't found the right solution)
    end
end

createCorrGraph(i, T1, T2, N1, N2, numStim)
end

