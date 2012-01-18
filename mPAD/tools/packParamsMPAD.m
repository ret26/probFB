function Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary);
  
  % function Params = PackParamsNew(G,Lam,Var,NonLin);
  %
  % Packs all the parameters into a single structure
  %
  % INPUTS
  % G,Lam1,Len2,Var1,Var2,Mu2,vary  
  %
  % OUTPUTS
  % Params = structure including
  %  G,Lam1,Len2,Var1,Var2,Mu2,vary 
  %
  % see test_packParamsMPAD.m for tests

Params.G = G;
Params.Lam1 = Lam1;
Params.Var1 = Var1;
Params.Len2 = Len2;
Params.Var2 = Var2;
Params.Mu2 = Mu2;
Params.vary = vary;