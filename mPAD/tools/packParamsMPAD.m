function Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary,varargin);
  
  % function Params = PackParamsNew(G,Lam1,Var1,Len2,Var2,Mu2,vary);
  % Packs all the parameters into a single structure - for AR(2)
  % carriers
  %
  % function Params = PackParamsNew(G,Lam1,Var1,Len2,Var2,Mu2,vary,om);
  % Packs all the parameters into a single structure - for pFB
  % carriers
  %
  % INPUTS
  % [G,Lam1,Var1,Len2,Var2,Mu2,vary] and optionally [om] 
  %
  % OUTPUTS
  % Params = structure including
  %  [G,Lam1,Len2,Var1,Var2,Mu2,vary] and optionally [om]  
  %
  % see test_packParamsMPAD.m for tests

Params.G = G;
Params.Lam1 = Lam1;
Params.Var1 = Var1;
Params.Len2 = Len2;
Params.Var2 = Var2;
Params.Mu2 = Mu2;
Params.vary = vary;

if nargin>7
  Params.om = varargin{1};
end