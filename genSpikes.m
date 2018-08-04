function [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT)
% Generates a binary spike timecourse to serve as "stimulus"
% Usage: [spkInds,spkVec] = genSpikes(stimDur,spkRate,dT)
%        [spkInds,spkVec] = genSpikes(1,30,0.001)

% Generate a binary spike train over stimulus period
spkVec = zeros(stimDur/dT + 1,1);
spkInts = floor(1/(dT*spkRate));

% Deterministic input spike train
spkInds = 1:spkInts:spkInts*spkRate;
spkVec(spkInds) = 1;

% Poisson input spike train?

end