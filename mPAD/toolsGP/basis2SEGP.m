function y = basis2SEGP(z,lenx,varx,tol)

% function y = basis2SEGP(z,lenx,varx,tol)
%
% converts the basis function representation of a squared
% exponential Gaussian Process to the GP itself
%
% INPUTS
% z = basis function coefficients [T,1]
% lenx = horizontal length-scale parameter of the GP
% varx = vertical scale parameter of GP (marginal variance)
% tol = tolerence that sets size of the basis functions (probably
%       should be > 9)
%
% OUTPUTS
% y = Gaussian process

tau = ceil(lenx/sqrt(2)*tol);
bh = sqrt(varx*sqrt(2)/(lenx*sqrt(pi))); % basis function height to
                                         % get correct variance
bas = bh*exp(-1/(lenx^2)*([-tau:1:tau]').^2);
  
y =  conv(z,bas,'same');
