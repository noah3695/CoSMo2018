function [spksOut] = LIFModel(spkInds,spkVec,gShared,gIndep,dT,plotOn)
% Leaky integrate and Fire Model
% Usage: [spksOut] = LIFModel(spkInds,spkVec,gShared,gIndep,dT,plotOn)
%        [spksOut] = LIFModel([spkTimes],[0101...],3,[],0.001,1)           
%
% - spksOut is a list of timestamps of output spikes
% - spkInds provides a vector of spike times
% - spkVec is a binary vector of length t/dt + 1
% - gShared is shared noise input
% - gIndep is independent noise that should be modulated by tuning
% - plotOn sets whether or not to plot neuron simulation
%
%   Note: Must generate spkInds,spkVec,gIndep with genSpikes.m and 
%         indNoise.m prior to running

% TO DO:
% - figure out noise distributions (balance of gains)
% - figure out correct parameters for gains of noise vs signal (biophys)

t = [0:dT:1]';

% Biophys Params
tau = 0.015;
R = 400000; % ohms
Iepsc = 0.00008; % amps
thresh = -20;   % Spike threshold voltage
vRMP = -65;     % resting membrane potential

% Iterate over time period
V(1) = vRMP;
V_AHP = -80;

numout = 0;
spksOut = 0;

for i = 2:numel(t)
    % Reset if previous timepoint spiked
    if V(i-1,1) >= -20
        V(i,1) = V_AHP;
        continue
    end
    
    % current @t = current @t-1 - leak_channel + spike input*transfer fxn +
    % normal noise + more normal noise
    V(i,1) = V(i-1,1) + (1/tau)*(vRMP - V(i-1,1))*dT + R*Iepsc*spkVec(i-1) + ...
        gShared*randn(1) + gIndep(i-1); % + comNoise() + indNoise();
    
    % Spike if threshold is reached
    if V(i,1) >= thresh
        V(i,1) = 40;
        
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

end