function [gIndep] = indNoise(gainNoise,dT,stimDur)
% Create a vector of voltage noise to serve as independent inputs to cells
% Usage: [gIndep] = indNoise(gainNoise,dT,stimDur)
%        [gIndep] = indNoise(3,0.001,1)

gIndep = gainNoise*randn(stimDur/dT+1,1);

end