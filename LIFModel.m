function [spksOut] = LIFModel(stimRate,gShared,gIndep,plotOn)
% Leaky integrate and Fire Model
% Usage: [spksOut] = LIFModel(stimRate,gShared,gIndep)
%
% - spksOut is a list of timestamps of output spikes
% - stimRate is input spike rate in sp/s
% - gShared is shared noise input
% - gIndep is independent noise that should be modulated by tuning
% - plotOn sets whether or not to plot neuron simulation

% TO DO:
% - figure out noise distributions
% - figure out correct parameters for gains of noise vs signal (biophys)

dT = 0.001;
t = 0:dT:1;

% Biophys Params
tau = 0.015;
R = 400000; % ohms
Iepsc = 0.00008; % amps
thresh = -20;   % Spike threshold voltage
vRMP = -65;     % resting membrane potential

% Stim Params
stimDur = 1;
[spkInds,spkVec] = genSpikes(stimDur,stimRate,dT);

% Iterate over time period
V(1) = vRMP;
V_AHP = -80;

numout = 0;
spksOut = 0;

for i = 2:numel(t)
    % Reset if previous timepoint spiked
    if V(i-1) >= -20
        V(i) = V_AHP;
        continue
    end
    
    % current @t = current @t-1 - leak_channel + spike input*transfer fxn +
    % normal noise + more normal noise
    V(i) = V(i-1) + (1/tau)*(vRMP - V(i-1))*dT + R*Iepsc*spkVec(i-1) + ...
        gShared*randn(1) + gIndep*randn(1); % + comNoise() + indNoise();
    
    % Spike if threshold is reached
    if V(i) >= thresh
        V(i) = 40;
        
        spksOut(numout+1) = i;
        numout = numout + 1;
    end
end

% Plot intracellular voltage timecourse and in/out spikes

if plotOn
    figure;
    hold on;
    plot(t,V);
    set(gca,'YLim',[-100 60],'XLim',[0 0.50]);
    scatter(t(spkInds),50*spkVec(spkInds),'.k');
    scatter(t(spksOut),45*ones(numel(spksOut),1),'.r');
end

%%%%%%%%%%%%%%%%%%%%%%  Local Functions  %%%%%%%%%%%%%%%%%%%%%%%%

    function [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT)
        % Generate a binary spike train over stimulus period
        spkVec = zeros(stimDur/dT,1);
        spkInts = floor(1/(dT*spkRate));
        
        % Deterministic input spike train
        spkInds = 1:spkInts:spkInts*spkRate;
        spkVec(spkInds) = 1;
        
        % Poisson input spike train?
        
    end

    function [] = comNoise()
        % may or may not use, depending on form of noise injection
        % specified in other papers
    end

    function [] = indNoise()
        
    end
end