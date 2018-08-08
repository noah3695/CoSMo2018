function acc=nn_classifier(meanData,trialData,partVecTrial,condVecTrial,trialVec); 
    nPart = size(meanData,3);  % Partitions are the 3 dimension 
    nCond = size (meanData,1); % Number of conditons 
    part = [1:nPart]; 
    correct = [];
    for n=1:nPart 
        trainIndx = find(part~=n); 
        testIndx =  find(part==n); 
        Mu_hat = mean(meanData(:,:,trainIndx),3); % Calculate the training means 
        
        % Now classify per trial for the left out run
        % first reshape to have a trial x neuron for test
        for c=1:length(unique(condVecTrial))
            for t=1:length(unique(trialVec))
                x = trialData(partVecTrial==testIndx & condVecTrial==c & trialVec==t)';
                % do the classification
                dist=x*x'-2*Mu_hat*x'+sum(Mu_hat.^2,2);
                [~,k]=min(dist);  % Record the classification
                if c==k
                    correct=[correct;1];
                else
                    correct=[correct;0];
                end
            end
        end        
    end;
    % Calculate the % correct
    acc = sum(correct(:))/numel(correct(:)); 