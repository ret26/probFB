function [Freqs,Spec] = getSpecARTau(lams,varx,NumFreqs,RngFreqs)

  % function [Freqs,Spec] = getSpecARTau(lams,varx,NumFreqs,RngFreqs)
  %
  % Computes the spectra of a single AR(tau) process
  %
  % INPUTS
  % lams = AR(tau) dynamical parameter [tau,1]
  % varx = AR(tau) dynamical noise 
  % NumFreqs = number of frequencies over which to return
  %            the spectra
  % RngFreqs = range over which to return the spectra
  %
  % OUTPUTS
  % Freqs = frequencies at which the spectrum is
  %         evaluated [1,NumFreqs]  
  % Spec = spectrum [1,NumFreqs]  
  
Freqs = [RngFreqs(1):diff(RngFreqs)/(NumFreqs-1):RngFreqs(2)];

Omegas = 2*pi*Freqs;

tau = length(lams);

%% OLD VERSION WITH NUMERICAL ISSUES
% c1 = 0;
% for t=1:tau
%   c1 = c1-2*lams(t)*exp(i*Omegas*t);
% end

% Taus = repmat([1:tau],[tau,1]);

% dtau = Taus-Taus';
% dtau = dtau(:);
% lamlam = lams'*lams;
% lamlam = lamlam(:);

% % Old
% % c2=0;

% % for t=1:length(lamlam)
% %   c2 = c2+lamlam(t)*cos(Omegas*dtau(t));
% % end

% % New is twice as fast 
% c2 = lams*lams';
% for t=1:tau
%   tt = [t+1:tau];
%   c2 = c2+2*lams(t)*lams(tt)*exp(i*(tt-t)'*Omegas);
% end

%Spec = varx./real(1+c1+c2+1e-10); 
%Spec = varx./real(1+c1+c2); 

% NEW VERSION THAT GUARANTEES POSITIVITY 
c3 = 1;
c4 = 0;
for t=1:tau
  c3 = c3-lams(t)*cos(Omegas*t);
  c4 = c4+lams(t)*sin(Omegas*t);
end

Spec = varx./(c3.^2+c4.^2); 
