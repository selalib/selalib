function [dz, ddz] = plasma_dispersion_function_derivative(zeta)

z = sqrt(pi)*exp(-zeta.^2).*(1i - erfi((zeta)));

dz = -2*(1+z.*zeta);

ddz = 4.*zeta + (4*zeta.^2-2).*z;


% function [z, dz, ddz] = plasma_dispersion_bot(xi)
% 
% if ( imag(xi) > 0 )
%	z = sqrt(pi)*exp(-xi.^2).*(-1i - erfi(xi));
%
%	dz = -2*(1-z.*xi);
% 
%	ddz = -4.*xi + (4*xi.^2+2).*z;
% else
%	z = sqrt(pi)*exp(-xi.^2).*(1i - erfi(xi));
%
%	dz = -2*(1+z.*xi);
% 
%	ddz = 4.*xi + (4*xi.^2-2).*z;
% end

