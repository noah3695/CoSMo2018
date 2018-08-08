%% Load a set of test data
numRunRep=10;
T.partVec = repmat([1:5]',[numNeuron*numRunRep*numStim,1]);

for i=1:10 %10 neurons
    N{i} = getrow(T,T.neuron==i);
    nNeuron(:,i) = N{i}.spikeNum; %concatonate to 1 big matrix
    nCorrect(:,i) = N{i}.prefDir; %correct matrix
    nLabel(:,i) =  N{i}.neuron; %the neuron label
end

%% Get the indices where 1 is right, 2 is right, etc
u = unique(nCorrect);
for i=1:length(u)
    indices{i} = find(nCorrect == u(i));
end

%but to start with for either this function or classify
%we'll need a 5x100x5 matrix (5 cond, 100 neurons, 5 runs).

%% leave one run out - assuming we have 5 runs X 10 trials per run
for p = 1:length(unique(T.partVec))
    
    train = getrow(T,T.partVec~=p);
    test  = getrow(T,T.partVec==p);
    
    leaveout = 1:5;
    leaveout(p) = [];
    trainList = getrow(T,T.partVec == leaveout(1));
    for iList=2:length(leaveout)
        trainList = cell2struct(cellfun(@vertcat,struct2cell(trainList),struct2cell(getrow(T,T.partVec == leaveout(iList))),'uni',0),fieldnames(trainList),1);
    end
    meanTrain = tapply(train,{'neuron','stimDir'},{'spikeNum','mean'});
    
    % apply the classifier - test set
    
    indicesTrain = [];
    indicesTest= [] ;
    indices = [];
    for i=1:length(u)
        foo = 1:numStim;
        foo(i) = []; %get rid of the training data
        indices{i} = find(nCorrect == u(i));
        indicesTrain{i} = indices{i}(i:numStim:end);
        indicesTest{i} = indices{i}(2:2:end);
    end
    
    
    
    
    %%% Find how probable is decoding using naive test? - use the gaussian
    %%% distribution / assume each neuron & time bin is independent
    
    %% copied from Konrad's tutorial
    % devianceTrain = [];
    % deviance = [];
    % for t = 1:nTrials %for all trials
    %     for j = 1:length(u) %for all the labels it could possibly be
    %         spikes = sp_count(:,:,t);
    %         %subtract from the spikes the means that is associated with the
    %         %hypothesis from that specific label
    %         deviance(t,j) = sum(sum((spikes-means(:,:,j)).^2 ./(2*stds(:,:,j).^2)));
    %     end
    % end
    %
    % %% we need to normalize
    % % for each of the trials, subtract the mean
    %
    % normDev = [];
    % for i = 1:length(deviance)
    %     foo = mean(deviance(i,:));
    %     normDev(i,:) = deviance(i,:) - foo;
    % end
    
    
end

