function [D,K,T] = unpackDimsMPAD(Dims)
  
  %  function [D,K,T] = unpackDimsMPAD(Dims)
  %
  % Unpack the dimensionality structure
  %
  % INPUT
  % Dims = structure containing
  % D,K,T
  %
  % OUTPUTS
  % D,K,T
  %
  % see test_packSimsMPAD.m for tests
    
D = Dims.D;
K = Dims.K;
T = Dims.T;
    