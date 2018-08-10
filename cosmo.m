function varargout = cosmo(what,varargin)
%hello world
% Define common stuff here
baseDir = fileparts(which('cosmo.m'));
dataDir = fullfile(baseDir,'data');

numPrefs    = 5;
stimLow     = 1;    % Define stimulus range
stimHigh    = 5;
dStim       = 0.25; % define stimulus steps
stims       = stimLow:dStim:stimHigh;
numStim     = numel(stims);
numNeuron   = 1000;    % number of neurons altogether
numRun      = 8;  % 8 runs
numRep      = 10; % 10 repetitions per run, (80) overall

popType = 'negative';

switch what
    case 'GEN_tunedPopulation'
        % generate tuning population of neurons
        % usage: cosmo('GEN_population','numNeuron',1000,'numStim',3);
        plotFig   = 1;
        vararginoptions(varargin,{'numNeuron','numPrefs','plotFig','scale','offset','sigma'});
        
        % determine preferred tuning and variance per neuron
        prefDir = randi(5,[numNeuron,1]);
        sigma   = 7*rand(numNeuron,1)+0.5;     % 0.4 former
        scale   = rand([numNeuron,1]);           % 1 former
        offset  = rand([numNeuron,1])*0.75 + 0.5;           % 0.5 former
        
        % organised preferred direction - pref: 1,2,...,numPrefs
        TC = @(scale,stim,prefDir,sigma,offset)...
                scale .* exp(-((stim-prefDir).^2)./sigma) + offset;
        
        for i = 1:numNeuron
            tuning(i,:) = TC(scale(i),stims,prefDir(i),sigma(i),offset(i));
        end
        
        rescale = max(max(tuning));
        tuning = tuning./rescale;
        scale = scale./rescale;
        offset = offset./rescale;
        
        if plotFig==1 % optional plotting of tuning functions across neurons
            figure;
            hold on;
            for i=1:numPrefs
                inds = prefDir==i;
                idx = find(inds==1);
                subplot(1,numPrefs,i)
                plot(stims,TC(scale(inds)),stims,prefDir(inds),sigma(inds),offset(inds));
               %plot(stims,TC(scale(idx(1:10)),stims,prefDir(idx(1:10)),sigma(idx(1:10)),offset(idx(1:10))),'LineWidth',2);
            end
        end
        % save the tuning matrix (numNeuron x numStim)
        save(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)),...
            'tuning','prefDir','sigma','scale','offset');
    
    case 'GEN_LIF'
        % define default parameters for LIFModel
        gShared     = 0.003; % shared noise
        gIndep      = 0.008; % independent noise
        gAnat       = 0.001;
        plotOn      = 0;
        stimDur     = 2; % in seconds
        dT          = 0.001; % time increment
        numRun      = 8;  % 8 runs
        numRep      = 10; % 10 repetitions per run, (80) overall
        spikeScale  = 80; % in Hz
        
        TC = @(scale,stim,prefDir,sigma,offset)...
                scale * exp(-((stim-prefDir).^2)./sigma) + offset; 
            
        vararginoptions(varargin,{'stimRate','gShared','gIndep','plotOn','numNeuron','numStim','dt','stimDur','numRep','spikeScale','popType'});
        
        TT=[]; % initialise for storage (spikes across neurons / stimuli)
        % load the correct tuning matrix
        D = load(fullfile(dataDir,sprintf('tunMatrix_%dneurons_%dstim',numNeuron,numStim)));
        switch popType
            case 'mixture'
                anatVec = repmat([1; -1],numNeuron,1);
            case 'positive'
                anatVec = repmat([1; 1],numNeuron,1);
            case 'negative'
                anatVec = repmat([-1; -1],numNeuron,1);
        end    
        for t=1:numStim
            for r=1:numRun
                for rep=1:numRep
                    gSharedVec = sharedNoise(gShared,dT,stimDur); % same across neurons
                    for n=1:numNeuron
                        % 1) determine spike rate based on tuning (pull
                        % from params
                        resp = TC(D.scale(n),stims(t),D.prefDir(n),D.sigma(n),D.offset(n));
%                         resp = D.tuning(n,t);
                        spkRate = resp*spikeScale;
                        
                        % 2) generate spikes
                        [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT);
                        % 3) add anatOff
                        anatSign = anatVec(n);
                        % 4) run the LIFModel
                        %T.spikes{1} = LIFModel(spkInds,spkVec,gSharedVec,0,gAnat,anatSign,dT,stimDur,plotOn);
                        T.spikes{1} = LIFModel(spkInds,spkVec,gSharedVec,gIndep,gAnat,anatSign,dT,stimDur,plotOn);
                        T.spikeNum  = numel(T.spikes{1});
                        T.neuron    = n;
                        T.numRun    = r;
                        T.numRep    = rep;
                        T.prefDir   = D.prefDir(n);
                        T.stimDir   = stims(t);
                        T.anatSign  = anatSign;
                        TT          = addstruct(TT,T);
                        clear spkVec spkRate spkInds;
                    end
                end
                fprintf('Generated all neurons for stimulus: %d/%d runs %d/%d\n',t,numStim,r,numRun);
            end
        end
        figure;
        lineplot(TT.stimDir,TT.spikeNum,'split',TT.prefDir,'style_thickline');
        xlabel('Direction'); ylabel('Spike number'); title('Responses split by preferred direction');
        % save the outputs
        save(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron,numStim,popType)),'-struct','TT');
    
    case 'PLOT_population'
        % plot the population

        vararginoptions(varargin,{'numNeuron','numStim','popType'});
        
        T = load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron,numStim,popType)));
        
        for l=1:numStim
            legLab{l} = sprintf('stim-%d',l);
        end
        figure
%         plt.hist(T.spikeNum,'split',T.prefDir);
        ylabel('number of spikes');
        
        % extract variance and mean
        T1=tapply(T,{'neuron','prefDir','stimDir'},{'spikeNum','mean','name','spikeNum_mean'},...
            {'spikeNum','var','name','spikeNum_var'});
        % plot responses in dependence of preferred - presented stimulus
        figure(1)
        subplot(1,2,1)
        %lineplot(abs(T1.prefDir-T1.stimDir),T1.spikeNum_mean,'split',T1.prefDir,...
        %        'leg',legLab,'style_thickline','markertype','o','markersize',12);
        plt.line(abs(T1.prefDir-T1.stimDir),T1.spikeNum_mean,'split',T1.prefDir,...
                'leg',legLab,'markertype','o','markersize',12);
      
        xlabel('Preferred-presented stimulus');
        ylabel('Mean spike number response');
        subplot(1,2,2)
        % calculate the variance
      %  lineplot(abs(T1.prefDir-T1.stimDir),T1.spikeNum_var,'split',T1.prefDir,...
       %     'leg',legLab,'style_thickline','markertype','o','markersize',12);
%         plt.line(abs(T1.prefDir-T1.stimDir),T1.spikeNum_var,'split',T1.prefDir,...
%             'leg',legLab,'markertype','o','markersize',12);
        
        ylabel('Variance in spike number across trials');
        
        % plot stuff
        figure(2)
        barplot(abs(T1.prefDir-T1.stimDir),T1.spikeNum_var,'split',T1.prefDir);
    
    case 'CALC_corr_dprime'

        vararginoptions(varargin,{'numNeuron','numStim','popType'});
        
        T = load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron,numStim,popType)));
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
                D.corrVar       = corr(T2.spikeNum_var(T2.neuron==i),T2.spikeNum_var(T2.neuron==j));
                for c=1:numStim
                    t1=getrow(T,T.neuron==i & T.stimDir==c);
                    t2=getrow(T,T.neuron==j & T.stimDir==c);
                    D.dprime(1,c)=dprime(t1.spikeNum,t2.spikeNum);
                end
                D.neuron1       = i;
                D.neuron2       = j;
                D.sameNeuron    = double(i==j);
                pref1           = T1.prefDir(T1.neuron==i);
                pref2           = T1.prefDir(T1.neuron==j);
                anat1           = T.anatSign(T.neuron==i);
                anat2           = T.anatSign(T.neuron==i);
                D.prefDir1      = pref1(1);
                D.prefDir2      = pref2(1);
                D.prefSame      = double(D.prefDir1==D.prefDir2);
                D.anatSign1     = anat1(1);
                D.anatSign2     = anat2(1);
                DD = addstruct(DD,D);
                Rm(i,j)=D.corrMean;
                Rv(i,j)=D.corrVar;
            end
            fprintf('Calc corr pairs:\tneuron %d/%d\n',i,numNeuron);
        end    
        save(fullfile(dataDir,sprintf('corr_neuronPairs_%dneurons_%sPopulation',numNeuron,popType)),'-struct','DD'); 
    
    case 'CALC_classify'   
        nNeuron=48;
        
        vararginoptions(varargin,{'numNeuron','numStim','popType'});        
        T = load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron,numStim,popType)));
        
        % add partVec
        T1 = tapply(T,{'stimDir','numRun','neuron'},{'spikeNum','mean'});        
        % rearrange
        for i=1:length(unique(T1.numRun))
            tmp = getrow(T1,T1.numRun==i);
            [indx,j,k] = pivottable([tmp.stimDir],[tmp.neuron],[tmp.spikeNum],'mean');
            data(:,:,i) = indx;
        end
        
        % choose subsets of neurons
        switch popType
            case 'mixture'
                % positive / negative / mixed corr
                indxN(:,1)=randsample(unique(T.neuron(T.anatSign==1)),nNeuron);
                indxN(:,2)=randsample(unique(T.neuron(T.anatSign==-1)),nNeuron);
                indxN(:,3)=[randsample(indxN(:,1),nNeuron/2); randsample(indxN(:,2),nNeuron/2)]';
            case 'positive'
                indxN(:,1)=randsample(unique(T.neuron),nNeuron)';
            case 'negative'
                indxN(:,1)=randsample(unique(T.neuron),nNeuron)';
        end
        
        % submit to classifier, distance calculation
        for s=1:size(indxN,2)
            for i=1:2 % give all stimuli or only 5
                if i==1 % then give all stimuli
                    data_sub = data(:,indxN(:,s),:);
                    T_sub    = getrow(T,ismember(T.neuron,indxN(:,s)));
                    acc(s,i) = nn_classifier(data_sub,T_sub.spikeNum,T_sub.numRun,T_sub.stimDir,T_sub.numRep);
                else
                    indxS    = ismember(stims,[1:numStim]);
                    data_sub = data(indxS,indxN(:,s),:);
                    T_sub    = getrow(T,ismember(T.stimDir,T.prefDir) & ismember(T.neuron,indxN(:,s)));
                    acc(s,i) = nn_classifier(data_sub,T_sub.spikeNum,T_sub.numRun,T_sub.stimDir,T_sub.numRep);
                end
                
            end
        end
        %         dist = rsa_distanceLDC(T1.spikeNum,T1.numRun,T1.stimDir);
        %         % plot distances
        %         figure
        %         imagesc(rsa_squareRDM(dist));
        keyboard;
        % save classifier results
        
    case 'SAVE_subPop'
        vararginoptions(varargin,{'popType','numNeuron'});
        
        T = load(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron,numStim,popType)));
        
        TT = [];
        newNeuron=[];
        allNeuron = 1:numNeuron;
        
        if strcmp(popType,'mixture') % consider also weights
            TT = getrow(T,ismember(T.neuron,[1:numNeuron/2]));
        else % only prefDir
            % newNeuron - save so that there is the right number
            for i=1:length(unique(T.prefDir))
                origN = unique(T.neuron(T.prefDir==i));
                randN = randsample(origN,floor(length(origN)/2));
                newNeuron = [newNeuron randN'];
                T2 = getrow(T,ismember(T.neuron,randN));
                TT=addstruct(TT,T2);
            end
            missingNew = numNeuron/2 - length(unique(TT.neuron));
            if missingNew > 0
                addN = randsample(allNeuron(~ismember(allNeuron,newNeuron)),missingNew);
                Tadd = getrow(T,ismember(T.neuron,addN));
                TT=addstruct(TT,Tadd);
            end
        end
        
        save(fullfile(dataDir,sprintf('LIF_%dneurons_%dstim_%sPopulation',numNeuron/2,numStim,popType)),'-struct','TT');
        
    case 'PLOT_corr'
%         popType='mixture';
        vararginoptions(varargin,{'numNeuron','popType'});
        T = load(fullfile(dataDir,sprintf('corr_neuronPairs_%dneurons_%sPopulation',numNeuron,popType)));
      
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
        title(sprintf('Signal correlation across neuron pairs - %s population',popType));
        subplot(1,3,2)
        imagesc(M_var);
        title('Noise correlation across neuron pairs');
        subplot(1,3,3)
        imagesc(M_mean-M_var);
        title('Difference signal-noise correlation');
            
        figure
        subplot(2,2,1)
        plt.scatter(T.corrMean,T.corrVar,'subset',T.sameNeuron==0 & T.prefSame==0);
        xlabel('Signal correlation'); ylabel('Noise correlation');
        title(sprintf('Neurons with diff preferred direction - %s population',popType));
        subplot(2,2,2)
        plt.scatter(T.corrMean,T.corrVar,'subset',T.sameNeuron==0 & T.prefSame==1);
        xlabel('Signal correlation'); ylabel('Noise correlation');
        title('Neurons with same preferred direction');
        subplot(2,2,3)
        plt.hist(T.corrMean,'subset',T.sameNeuron==0);
        hold on;
        drawline(mean(T.corrMean(T.sameNeuron==0)),'dir','vert','color',[1 0 0]);
        title('Distribution of signal correlation');
        subplot(2,2,4)
        plt.hist(T.corrVar,'subset',T.sameNeuron==0);
        hold on;
        drawline(mean(T.corrVar(T.sameNeuron==0)),'dir','vert','color',[1 0 0]);
        title('Distribution of noise correlation'); 

    case 'CALC_fisherInfo'
%         popType = 'mixture';
        numRpts = numRun*numRep;
        
        [FI_corr,dprimo] = fisherInfo(dataDir,numNeuron,numRpts,numStim,stims,popType);
        
        save(fullfile(dataDir,sprintf('fisher_%dneurons_%dstim_%sPopulation',...
            numNeuron,numStim,popType)),'FI_corr','dprimo');
    otherwise
        fprintf('No such case\n');

end

% Local functions