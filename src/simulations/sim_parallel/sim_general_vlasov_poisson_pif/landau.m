function [ fun, J ] = landau( coeff, t )
%UNTITLED3 Summary of this function goes here
%   data(1) 

A=coeff(1);
omegat=coeff(2);
phi=coeff(3);
omega=coeff(4);

syms A t omegat phi omega


fun=A*log(cos(t.*omegat + phi).^2.*exp(-omega.*t));

if (nargout >1 )
%J=jacobian(fun,[A, omegat, phi, omega])
%matlabFunction(J)
J=[log(exp(-omega.*t).*cos(phi+omegat.*t).^2),...
    (A.*t.*sin(phi+omegat.*t).*-2.0)./cos(phi+omegat.*t),...
    (A.*sin(phi+omegat.*t).*-2.0)./cos(phi+omegat.*t),...
    -A.*t];
end



% fun=A*cos(t.*omegat + phi).^2.*exp(-omega.*t);
% 
% if (nargout >1 )
% %J=jacobian(fun,[A, omegat, phi, omega])
% J=[exp(-omega.*t).*cos(phi+omegat.*t).^2, ...
%      A.*t.*exp(-omega.*t).*cos(phi+omegat.*t).*sin(phi+omegat.*t).*-2.0,...
%      A.*exp(-omega.*t).*cos(phi+omegat.*t).*sin(phi+omegat.*t).*-2.0,...
%      -A.*t.*exp(-omega.*t).*cos(phi+omegat.*t).^2];
% end


end

