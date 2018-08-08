function d = dprime(respN1,respN2)
% function d = dprime(respN1,respN2)
% takes responses from neurons N1 and N2
% responses per condition across all trials

d = (mean(respN1)-mean(respN2))/sqrt(0.5*(std(respN1)+std(respN2)));

end