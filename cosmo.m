function varargout = cosmo(what,varargin)

% Define common stuff here
baseDir = '/Users/Eva/Documents/Data/CoSMo2018'; % differs per person
dataDir = fullfile(baseDir,'data');

switch what
    case 'GEN_tuningNeuron'
        % generate tuning population of neurons
        % usage: cosmo('GEN_population','numNeuron',1000,'numStim',3);
        numNeuron = 500; % tuned neurons PER(!) stimulus
        numStim   = 4;
        sigma     = 0.25;
        scale     = 2; % scaling function
        offset    = 0.5;
        plotFig   = 1;
        vararginoptions(varargin,{'numNeuron','numStim','plotFig','scale','offset','sigma'});
        
        % determine preferred tuning
        prefDir = randi(numStim,[numNeuron,1]);
        tuning = (scale*exp(-([1:numStim]-prefDir).^2)./(2*sigma^2))+offset;      
        if plotFig==1 % optional plotting of tuning function
            figure
            hold on;
            for i=1:numStim
                subplot(1,numStim,i)
                plot([1:numStim],tuning(prefDir==i,:));
            end
        end
        varargout{1}=tuning;
        save(fullfile(dataDir,sprintf('tuningFunc_%dstim',numStim)),'tuning');
    case 'GEN_LIF'
        % define inputs for the LIF
        % scaling parameters
        stimRate = 80;
        gShared  = 3;
        gIndep   = 2;
        plotOn   = 1;
        vararginoptions(varargin,{'stimRate','gShared','gIndep','plotOn'});
        % run the model
        D=LIFModel(stimRate,gShared,gIndep,plotOn);
        % save the outputs
    case 'PLOT_LIF'
    otherwise
        fprintf('No such case\n');

end

% Local functions