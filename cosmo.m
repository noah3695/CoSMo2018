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
        sigma     = 0.25;
        scale     = 2; % scaling function
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
        % define inputs for the LIF
        % scaling parameters
        stimRate = 80;
        gShared  = 3;
        gIndep   = 2;
        plotOn   = 0;
        stimDur  = 2; % in seconds
        dT = 0.001; % time increment
        numNeuron = 500;
        numStim = 10;
        spikeScale=30;
        vararginoptions(varargin,{'stimRate','gShared','gIndep','plotOn','numNeuron','numStim','dt','stimDur'});
        
        TT=[]; % initialise for storage (spikes across neurons / stimuli)
        % load the correct tuning matrix
        D = load(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)));
        for t=1:numStim
            for n=1:numNeuron
                % 1) determine spike rate based on tuning
                spkRate = D.tuning(n,t)*spikeScale;
                % 2) generate spikes
                [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT);
                % 3) run the LIFModel
                T.spikes = LIFModel(spkInds,spkVec,gShared*ones(size(spkVec)),gIndep,dT,plotOn);
                T.neuron = n;
                T.prefDir = D.prefDir(n);
                T.stimDir = t;
                TT = addstruct(TT,T);
            end
        end
        keyboard;
        % save the outputs
        save(fullfile('LIF_%dneurons_%dstim'),'-struct','TT');
    case 'PLOT_LIF'
    otherwise
        fprintf('No such case\n');

end

% Local functions