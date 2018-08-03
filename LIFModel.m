function [] = LIFModel()
% Leaky integrate and Fire Model
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
stimRate = 80;   % ultimately want to be pulled from tuning fxn
[spkInds,spkVec] = genSpikes(stimDur,stimRate,dT);

% Iterate over time period
V(1) = vRMP;
V_AHP = -80;

% spksOut = zeros(length(spkVec),1);
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
        3*randn(1) + 2*randn(1); % + comNoise() + indNoise();
    
    % Spike if threshold is reached
    if V(i) >= thresh
        V(i) = 40;
        
        spksOut(numout+1) = i;
        numout = numout + 1;
    end
end

figure;
hold on;
plot(t,V);
set(gca,'YLim',[-100 60],'XLim',[0 0.50]);
scatter(t(spkInds),50*spkVec(spkInds),'.k');
scatter(t(spksOut),45*ones(numel(spksOut),1),'.r');

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

function [] = epsp()

end

function [] = comNoise()

end

function [] = indNoise()

end
end