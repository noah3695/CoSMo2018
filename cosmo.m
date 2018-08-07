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
        prefDir = randi(numStim,[numNeuron,1]);
        sigma = abs(rand(numNeuron,1));
        gIndep = abs(rand(numNeuron,1));
        tuning = (scale*exp(-([1:numStim]-prefDir).^2)./(2*sigma.^2))+offset;      
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
        save(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)),'tuning','prefDir','gIndep');
    case 'GEN_LIF'
        % define default parameters for LIFModel
        gShared     = 0.05; % shared noise
        gIndep      = 0.001; % independent noise
        plotOn      = 0;
        stimDur     = 2; % in seconds
        dT          = 0.001; % time increment
        numNeuron   = 100;
        numStim     = 5;
        numRep      = 50;
        spikeScale  = 30; % in Hz
        vararginoptions(varargin,{'stimRate','gShared','gIndep','plotOn','numNeuron','numStim','dt','stimDur','numRep','spikeScale'});
        
        TT=[]; % initialise for storage (spikes across neurons / stimuli)
        % load the correct tuning matrix
        D = load(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)));
        gSharedVec = sharedNoise(gShared,dT,stimDur);
        for t=1:numStim
            for r=1:numRep
                Sign=(randi([0,1],numNeuron,1)*2-1)*0.05; %positive or negative
                for n=1:numNeuron
                    % 1) determine spike rate based on tuning
                    spkRate = D.tuning(n,t)*spikeScale;
                    % 2) generate spikes
                    [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT);
                    % 3) run the LIFModel
                    T.spikes{1} = LIFModel(spkInds,spkVec,gSharedVec,D.gIndep(numNeuron)*Sign(numNeuron),dT,stimDur,plotOn);
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
        % choose subsets of generated population
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
        for i=1:numNeuron
            for j=i:numNeuron
                D.corrMean  = corr(T1.spikeNum_mean(T1.neuron==i),T1.spikeNum_mean(T1.neuron==j));
                %D.corrVar   = corr(T1.spikeNum_var(T1.neuron==i),T1.spikeNum_var(T1.neuron==j));
                D.corrVar   = corr(T.spikeNum(T.neuron==i),T.spikeNum(T.neuron==j));
                D.neuron1   = i;
                D.neuron2   = j;
                pref1 = T1.prefDir(T1.neuron==i);
                pref2 = T1.prefDir(T1.neuron==j);
                D.prefDir1  = pref1(1);
                D.prefDir2  = pref2(1);
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
        
        T.prefSame = zeros(size(T.neuron1));
        T.prefSame(T.prefDir1==T.prefDir2)=1;
        T.sameNeuron = zeros(size(T.neuron1));
        T.sameNeuron(T.neuron1==T.neuron2)=1;
      
        figure
        scatterplot(T.corrMean,T.corrVar,'split',T.prefSame,'leg',{'tuning same','tuning different'},'subset',T.sameNeuron==0);
        %plt.scatter(T.corrMean,T.corrVar,'split',T.prefSame,'leg',{'tuning same','tuning different'});
        xlabel('Signal correlation');
        ylabel('Noise correlation');
        title('Correlation structure across neuron pairs');
        
        % make into matrix
        M_mean  = rsa_squareIPM(T.corrMean');
        M_var   = rsa_squareIPM(T.corrVar');
        
        for f=1:2 % do for same / diff prefDirection
            % choose neuron pairs with the same preferred direction
            if f==1
                lab = 'different';
            else
                lab = 'same';
            end
            indx1 = T.neuron1(find(T.prefSame==f-1 & T.sameNeuron==0));
            indx2 = T.neuron2(find(T.prefSame==f-1 & T.sameNeuron==0));
            
            % make smaller matrices
            M1_mean = M_mean(indx1,indx2);
            M1_var  = M_var(indx1,indx2);
            figure
            subplot(2,3,1)
            imagesc(M1_mean)
            title(sprintf('Mean corr structure %s preference tuning',lab));
            subplot(2,3,2)
            imagesc(M1_var)
            title('Var corr structure');
            subplot(2,3,3)
            imagesc(M1_mean-M1_var)
            title('Difference');
            subplot(2,3,4)
            scatterplot(M1_mean(:),M1_var(:));
            hold on
            drawline(0,'dir','horz');
            drawline(0,'dir','vert');
            xlabel('Mean corr'); ylabel('Noise corr');
            subplot(2,3,5)
            histogram(M1_mean(:));
            title('Mean corr across neuron pairs');
            subplot(2,3,6)
            histogram(M1_var(:));
            title('Noise corr across neuron pairs');
            clear indx1 indx2 M1_mean M1_var
        end
    otherwise
        fprintf('No such case\n');

end

% Local functions