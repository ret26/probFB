function [G,Lam1,Var1,Len2,Var2,Mu2,vary] = unpackParamsMPAD(Params);
  
  % [G,Lam1,Var1,Len2,Var2,Mu2,vary] = unpackParamsMPAD(Params);
  %
  % Unpacks all of the parameters out of the the structure
  %
  % INPUTS
  % Params = structure including
  %  G,Lam1,Len2,Var1,Var2,Mu2,vary 
  %
  % OUTPUTS
  % G,Lam1,Len2,Var1,Var2,Mu2,vary
  %
  % see test_packParamsMPAD.m for tests

G = Params.G;
Lam1 = Params.Lam1;
Var1 = Params.Var1;
Len2 = Params.Len2;
Var2 = Params.Var2;
Mu2 = Params.Mu2;
vary = Params.vary;

    
    
    