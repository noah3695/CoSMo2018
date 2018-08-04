function [spksOut,V] = LIFModel(spkInds,spkVec,gShared,gIndep,dT,stimDur,plotOn)
% Leaky integrate and Fire Model
% Usage: [spksOut] = LIFModel(spkInds,spkVec,gShared,gIndep,dT,plotOn)
%        [spksOut] = LIFModel([spkTimes],[0101...],3,[],0.001,1)
%
% - spksOut is a list of timestamps of output spikes
% - spkInds provides a vector of spike times
% - spkVec is a binary vector of length t/dt + 1
% - gShared is shared noise input (a vector for all cells in pop)
% - gIndep is independent noise that should be modulated by tuning
% - plotOn sets whether or not to plot neuron simulation
%
%   Note: Must generate spkInds,spkVec,gIndep with genSpikes.m and
%         indNoise.m prior to running

% TO DO:
% - figure out noise distributions (balance of gains)
% - figure out correct parameters for gains of noise vs signal (biophys)

t = [0:dT:stimDur]';

% Biophys Params (Wong & Wang 2006)
tau = 0.020; % ---- NS changed (original 0.015)
Vthresh = -0.050;   % Spike threshold voltage (original -20)
vRMP = -0.070;     % resting membrane potential

% Iterate over time period
V(1) = vRMP;
% V_AHP = -0.055; % from wong & wang + 2ms lockout period)
V_AHP = -0.080; % original -80

numout = 0;
spksOut = 0;

for i = 2:numel(t)
    % Reset if previous timepoint spiked
    if V(i-1,1) >= Vthresh
        V(i,1) = V_AHP;
        
        continue
    end
    
    % current @t = current @t-1 - leak_channel + spike input*transfer fxn +
    % independent noise + shared noise
    V(i,1) = V(i-1,1) - (1/tau)*(V(i-1,1) - vRMP)*dT -0.006*(1/tau)*V(i-1,1)*spkVec(i-1) + ...
        gIndep*randn(1) + gShared(i-1);
    
    % Spike if threshold is reached
    if V(i,1) >= Vthresh
        V(i,1) = 0.040;
        
        spksOut(numout+1) = i;
        numout = numout + 1;
    end
end

% Plot intracellular voltage timecourse and in/out spikes

if plotOn
    figure;
    hold on;
    plot(t,V*1000);
    ylabel('Voltage (mV)');
    xlabel('Time (sec)');
    set(gca,'YLim',[-100 60],'XLim',[0 1]);
    scatter(t(spkInds),50*spkVec(spkInds),'.k');
    if spksOut ~= 0
        scatter(t(spksOut),45*ones(numel(spksOut),1),'.r');
    end
end

end