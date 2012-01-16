function [Freqs,Spec,fMAX,SpecMAX,dF1,dF2] = getAR2Spec(lams,varx,NumFreqs,RngFreqs)

  % function [Freqs,Spec,fMAX,SpecMAX,dF1,dF2] = ...
  %           getAR2Spec(lams,varx,NumFreqs,RngFreqs)
  %
  % Computes the spectra of a single AR(2) process
  %
  % INPUTS
  % lams = AR(2) dynamical parameter [2,1]
  % varx = AR(2) dynamical noise 
  % NumFreqs = number of frequencies over which to return
  %            the spectra
  % RngFreqs = range over which to return the spectra
  %
  % OUTPUTS
  % Freqs = frequencies at which the spectrum is
  %         evaluated [1,NumFreqs]  
  % Spec = spectrum [1,NumFreqs]  
  % fMAX = frequency at which maximum in spectrum occurs
  % SpecMax = maximum value of spectrum 
  % dF1 = Right half bandwidth
  % dF2 = Left half bandwidth  
  
Freqs = [RngFreqs(1):diff(RngFreqs)/(NumFreqs-1): ...
         RngFreqs(2)];
Omegas = 2*pi*Freqs;

c1 = 1 + lams*lams';
c2 = 2 * lams(1)*(lams(2)-1);
c3 = -2*lams(2);

Spec = varx./(c1+c2*cos(Omegas)+c3*cos(2*Omegas));

fMAX = real(acos(-c2/(4*c3))/(2*pi));
SpecMAX = varx./(c1+c2*cos(fMAX*2*pi)+c3*cos(2*fMAX*2*pi));

%temp = 1/(4*c3)*sqrt(c2^2-8*c3*(c2^2/(4*c3)-c1+c3));
temp = 1/4*sqrt(8*c1/c3 - c2^2/c3^2 - 8);

%temp1 = 1/4*sqrt(-8-4/lams(2)*(1+lams(1)^2+lams(2)^2)- ...
%                lams(1)^2/lams(2)^2*(lams(2)-1)^2);

dF1 = acos(-c2/(4*c3) + temp)/(2*pi);
dF2 = acos(-c2/(4*c3) - temp)/(2*pi);

% To deal with edge effects
if isreal(dF1)==0
  dF1 = 0;
end

if isreal(dF2)==0
  dF2 = 1/2;
end