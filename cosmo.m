function varargout = cosmo(what,varargin)

% Define common stuff here
baseDir = fileparts(which('cosmo.m'));
dataDir = fullfile(baseDir,'data');

switch what
    case 'GEN_tunedPopulation'
        % generate tuning population of neurons
        % usage: cosmo('GEN_population','numNeuron',1000,'numStim',3);
        numNeuron = 500; % number of neurons altogether
        numStim   = 4;
        sigma     = 0.4;
        scale     = 1; % scaling function
        offset    = 0.5;
        plotFig   = 1;
        vararginoptions(varargin,{'numNeuron','numStim','plotFig','scale','offset','sigma'});
        
        % determine preferred tuning
        prefDir = randi(numStim,[numNeuron,1]);
        tuning = (scale*exp(-([1:numStim]-prefDir).^2)./(2*sigma^2))+offset;      
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
        save(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)),'tuning','prefDir');
    case 'GEN_LIF'
        % define default parameters for LIFModel
        gShared     = 0.02; % shared noise
        gIndep      = 0.01; % independent noise
        plotOn      = 0;
        stimDur     = 2; % in seconds
        dT          = 0.001; % time increment
        numNeuron   = 500;
        numStim     = 10;
        numRep      = 50;
        spikeScale  = 30; % in Hz
        vararginoptions(varargin,{'stimRate','gShared','gIndep','plotOn','numNeuron','numStim','dt','stimDur','numRep'});
        
        TT=[]; % initialise for storage (spikes across neurons / stimuli)
        % load the correct tuning matrix
        D = load(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)));
        for t=1:numStim
            for r=1:numRep
                for n=1:numNeuron
                    % 1) determine spike rate based on tuning
                    spkRate = D.tuning(n,t)*spikeScale;
                    % 2) generate spikes
                    [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT);
                    % 3) run the LIFModel
                    T.spikes{1} = LIFModel(spkInds,spkVec,gShared*ones(size(spkVec)),gIndep,dT,plotOn);
                    T.spikeNum  = numel(T.spikes{1});
                    T.neuron    = n;
                    T.prefDir   = D.prefDir(n);
                    T.stimDir   = t;
                    TT          = addstruct(TT,T);
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
        numNeuron = 500;
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
        
        figure(3)
    case 'CHOOSE_subsets' 
        
    case 'PLOT_LIF'
    otherwise
        fprintf('No such case\n');

end

% Local functions