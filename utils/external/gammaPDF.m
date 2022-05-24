function y = gammaPDF(t,tau,n)
% Gamma function - used for impulse response calculations and HIRF
%
%   y = gammaPDF(t,tau,n)
%
% The gamma is essentially the convolution of n exponentials with time
% constant of tau. The peak of the function is near tau * n
%
% See various places, but Boynton & Heeger 1996, eq 3 is a good place.

if notDefined('t'), error('Time steps required'); end
if notDefined('tau'), tau = 1; end
if notDefined('n'), n = 2; end

% After Boynton et al.  See boyntonHIRF in mrVista
y = (t ./ tau).^(n-1) .* exp(-t ./ tau) / (tau*factorial(n - 1));
y = y/sum(y);

end
