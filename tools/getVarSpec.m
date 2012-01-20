function freqSq=getVarSpec(X)

% function freqSq=getVarSpec(X)
%
% Compute the second moment of a spectrum that is:
%
% if s(k) = normalised power spectrum of x
% freqSq = \sum_k s(k) freq(k)^2 
%
% INPUTS 
% X = signal to compute second moment of spectrum size [T,D]
%
% OUTPUTS
% freqSq = average square frequencies [1,D]
%
% for tests see test_getVarSpec.m

[T,D] = size(X); 

spec=abs(fft(X)).^2; % spectra
normSpec = spec.*(ones(T,1)*(1./sum(spec,1))); % normalised to unity

freq = linspace(0,1/2,floor(T/2)+1);
if mod(T,2)==0
  freq = [freq';-freq(end-1:-1:2)'];
else
  freq = [freq';-freq(end:-1:2)'];
end

freqSq = sum(normSpec.*(freq.^2*ones(1,D)),1);

%freqSq = 2*sum(normSpec(1:floor(T/2)+1,:).*(freq(1:floor(T/2)+1).^2*ones(1,D)),1);
