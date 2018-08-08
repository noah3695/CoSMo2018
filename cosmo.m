function varargout = cosmo(what,varargin)

% Define common stuff here
baseDir = fileparts(which('cosmo.m'));
dataDir = fullfile(baseDir,'data');

switch what
    case 'GEN_tunedPopulation'
        % generate tuning population of neurons
        % usage: cosmo('GEN_population','numNeuron',1000,'numStim',3);
        numNeuron = 100; % number of neurons altogether
        numStim   = 5;
        sigma     = 0.4;
        scale     = 1; % scaling function
        offset    = 0.5;
        plotFig   = 1;
        vararginoptions(varargin,{'numNeuron','numStim','plotFig','scale','offset','sigma'});
        
        % determine preferred tuning and variance per neuron
        prefDir = kron([1:numStim]',[ones(numNeuron/numStim,1)]);
        % organised preferred direction - pref: 1,2,...,numStim
        sigma   = abs(rand(numNeuron,1));
        %gIndep = abs(rand(numNeuron,1));
       % tuning = (scale*exp(-([1:numStim]-prefDir).^2)./(2*sigma.^2))+offset;
        tuning = (scale*exp(-([1:numStim]-prefDir).^2)./(2*sigma.^2))+offset;
        tuning = bsxfun(@rdivide,tuning,max(tuning,[],2));
        gIndep = tuning.*1.5;
        tuneScale = abs(randn(numNeuron,1));
        if plotFig==1 % optional plotting of tuning functions across neurons
            figure
            hold on;
            for i=1:numStim
                subplot(1,numStim,i)
                plot([1:numStim],tuning(prefDir==i,:));
            end
        end
        %varargout{1}=tuning;
        % save the tuning matrix (numNeuron x numStim)
        save(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)),'tuning','prefDir','gIndep','tuneScale');
    case 'GEN_LIF'
        % define default parameters for LIFModel
        gShared     = 0.05; % shared noise
        gIndep      = 0.001; % independent noise
        gAnat       = 0.03;
        plotOn      = 0;
        stimDur     = 2; % in seconds
        dT          = 0.001; % time increment
        numNeuron   = 100;
        sigmaAnat   = 15;
        numStim     = 5;
        numRep      = 50;
        spikeScale  = 30; % in Hz
        vararginoptions(varargin,{'stimRate','gShared','gIndep','plotOn','numNeuron','numStim','dt','stimDur','numRep','spikeScale'});
        
        TT=[]; % initialise for storage (spikes across neurons / stimuli)
        % load the correct tuning matrix
        D = load(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)));
        %anatVec = 1.5*exp(-([1:numNeuron]-round(numNeuron/2)).^2/(2*sigmaAnat.^2))+0.1;
        anatVec = repmat([1; -1], numNeuron, 1);
        
        for t=1:numStim
            for r=1:numRep
                gSharedVec = sharedNoise(gShared,dT,stimDur); % same across neurons
                Sign=(randi([0,1],numNeuron,1)*2-1)*0.5; %positive or negative
                for n=1:numNeuron
                    % 1) determine spike rate based on tuning
                    spkRate = D.tuning(n,t)*spikeScale*D.tuneScale(n);
                    % 2) generate spikes
                    [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT);
                    % 3) add anatOff 
                    anatOff = anatVec(n);
                    % 4) run the LIFModel
                    T.spikes{1} = LIFModel(spkInds,spkVec,gSharedVec,D.gIndep(n,t)*Sign(numNeuron),gAnat,anatOff,dT,stimDur,plotOn);
                    %T.spikes{1} = LIFModel(spkInds,spkVec,gSharedVec,gIndep,dT,stimDur,plotOn);
                    T.spikeNum  = numel(T.spikes{1});
                    T.neuron    = n;
                    T.prefDir   = D.prefDir(n);
                    T.stimDir   = t;
                    TT          = addstruct(TT,T);
                    clear spkVec spkRate spkInds;
                end
                fprintf('Generated all neurons for stimulus: %d/%d repetition %d/%d\n',t,numStim,r,numRep);  
            end
        end
        lineplot(TT.stimDir,TT.spikeNum,'split',TT.prefDir,'style_thickline');
        xlabel('Direction'); ylabel('Spike number'); title('Responses split by preferred direction');
        % save the outputs
        save(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim',numNeuron,numStim)),'-struct','TT');
    case 'PLOT_population'
        % plot the population
        numNeuron = 100;
        numStim = 5;
        vararginoptions(varargin,{'numNeuron','numStim'});
        
        T = load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim',numNeuron,numStim)));
        for l=1:numStim
            legLab{l} = sprintf('stim-%d',l);
        end
        
        % extract variance and mean
        T1=tapply(T,{'neuron','prefDir','stimDir'},{'spikeNum','mean','name','spikeNum_mean'},...
            {'spikeNum','var','name','spikeNum_var'});
        % plot responses in dependence of preferred - presented stimulus
        figure(1)
        subplot(1,2,1)
        lineplot(abs(T1.prefDir-T1.stimDir),T1.spikeNum_mean,'split',T1.prefDir,...
                'leg',legLab,'style_thickline','markertype','o','markersize',12);
        xlabel('Preferred-presented stimulus');
        ylabel('Mean spike number response');
        subplot(1,2,2)
        % calculate the variance
        lineplot(abs(T1.prefDir-T1.stimDir),T1.spikeNum_var,'split',T1.prefDir,...
                'leg',legLab,'style_thickline','markertype','o','markersize',12);
        ylabel('Variance in spike number across trials');
        
        % plot stuff
        figure(2)
        barplot(abs(T1.prefDir-T1.stimDir),T1.spikeNum_var,'split',T1.prefDir);
    case 'CALC_corr'
        numNeuron = 100;
        numStim = 5;
        vararginoptions(varargin,{'numNeuron','numStim'});
        
        T = load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim',numNeuron,numStim)));
        % extract variance and mean
        T1=tapply(T,{'neuron','prefDir','stimDir'},{'spikeNum','mean','name','spikeNum_mean'},...
            {'spikeNum','var','name','spikeNum_var'});
        DD=[];
        % subtract the condition mean across trials        
        NN=[];
        for sd=1:numStim
            for n=1:numNeuron
                N1=getrow(T,T.neuron==n&T.stimDir==sd);
                N2=getrow(T1,T1.neuron==n&T1.stimDir==sd);
                N1.spikeNum=N1.spikeNum-N2.spikeNum_mean;
                NN=addstruct(NN,N1);
            end
        end
        T2=tapply(NN,{'neuron','prefDir','stimDir'},{'spikeNum','mean','name','spikeNum_mean'},...
            {'spikeNum','var','name','spikeNum_var'});
        for i=1:numNeuron
            for j=i:numNeuron
                D.corrMean      = corr(T1.spikeNum_mean(T1.neuron==i),T1.spikeNum_mean(T1.neuron==j));
                %D.corrVar   = corr(T1.spikeNum_var(T1.neuron==i),T1.spikeNum_var(T1.neuron==j));
                %D.corrVar   = corr(NN.spikeNum(NN.neuron==i),NN.spikeNum(NN.neuron==j));
                D.corrVar       = corr(T2.spikeNum_var(T2.neuron==i),T2.spikeNum_var(T2.neuron==j));
                D.neuron1       = i;
                D.neuron2       = j;
                D.sameNeuron    = double(i==j);
                pref1           = T1.prefDir(T1.neuron==i);
                pref2           = T1.prefDir(T1.neuron==j);
                D.prefDir1      = pref1(1);
                D.prefDir2      = pref2(1);
                D.prefSame      = double(D.prefDir1==D.prefDir2);
                DD = addstruct(DD,D);
                Rm(i,j)=D.corrMean;
                Rv(i,j)=D.corrVar;
            end
            fprintf('Calc corr pairs:\tneuron %d/%d\n',i,numNeuron);
        end    
        save(fullfile(dataDir,sprintf('corr_neuronPairs_%dneurons',numNeuron)),'-struct','DD'); 
    case 'PLOT_corr'
        numNeuron=100;
        vararginoptions(varargin,{'numNeuron'});
        T = load(fullfile(dataDir,sprintf('corr_neuronPairs_%dneurons',numNeuron)));
      
        figure
        scatterplot(T.corrMean,T.corrVar,'split',T.prefSame,'leg',{'tuning different','tuning same'},'subset',T.sameNeuron==0);
        %plt.scatter(T.corrMean,T.corrVar,'split',T.prefSame,'leg',{'tuning same','tuning different'});
        xlabel('Signal correlation');
        ylabel('Noise correlation');
        title('Correlation structure across neuron pairs');
        
        % make into matrix
        M_mean  = rsa_squareIPM(T.corrMean');
        M_var   = rsa_squareIPM(T.corrVar');
        figure
        subplot(1,3,1)
        imagesc(M_mean);
        title('Signal correlation across neuron pairs');
        subplot(1,3,2)
        imagesc(M_var);
        title('Noise correlation across neuron pairs');
        subplot(1,3,3)
        imagesc(M_mean-M_var);
        title('Difference signal-noise correlation');
            
        figure
        subplot(2,2,1)
        scatterplot(T.corrMean,T.corrVar,'subset',T.sameNeuron==0 & T.prefSame==0);
        xlabel('Signal correlation'); ylabel('Noise correlation');
        title('Neurons with diff preferred direction');
        subplot(2,2,2)
        scatterplot(T.corrMean,T.corrVar,'subset',T.sameNeuron==0 & T.prefSame==1);
        xlabel('Signal correlation'); ylabel('Noise correlation');
        title('Neurons with same preferred direction');
        subplot(2,2,3)
        histplot(T.corrMean,'subset',T.sameNeuron==0);
        hold on;
        drawline(mean(T.corrMean(T.sameNeuron==0)),'dir','vert','color',[1 0 0]);
        title('Distribution of signal correlation');
        subplot(2,2,4)
        histplot(T.corrVar,'subset',T.sameNeuron==0);
        hold on;
        drawline(mean(T.corrVar(T.sameNeuron==0)),'dir','vert','color',[1 0 0]);
        title('Distribution of noise correlation'); 
      
    case 'CHOOSE_subset'
        numNeuron = 100;
        vararginoptions(varargin,{'numNeuron'});
        
        T = load(fullfile(dataDir,sprintf('corr_neuronPairs_%dneurons',numNeuron)));
        
        % randsample one neuron, then take the most anticorrelated, then
        % more anticorrelated neurons
        [sortMean,indxMean]=sort(T.corrMean);
        [sortVar,indxVar]=sort(T.corrVar);
        keyboard;
        
    otherwise
        fprintf('No such case\n');

end

% Local functions