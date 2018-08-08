function [SSR] = fitnormal(par,vals,pcount)
% Fit normal distribution to some histogram, report SSE

normdist = @(x,sigma,mu,scale) (scale)*exp(-((x-mu).^2)./(2*sigma^2));

y = normdist(vals,par(1),par(2),par(3));

residuals = pcount-y;

SSR = sum(residuals.^2);

end