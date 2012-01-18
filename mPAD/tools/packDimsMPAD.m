function Dims = packDimsMPAD(D,K,T);
  
  % function Dims = packDimsMPAD(D,K,T);
  %
  % Packs up all the dimensionality information into a
  % single structure
  %
  % INPUT
  % D,K,T
  %
  % OUTPUTS
  % Dims = structure containing
  % D,K,T
  %
  % see test_packSimsMPAD.m for tests


Dims.D = D;
Dims.K = K;
Dims.T = T;