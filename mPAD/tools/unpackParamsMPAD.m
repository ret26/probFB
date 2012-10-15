function [G,Lam1,Var1,Len2,Var2,Mu2,vary,varargout] = unpackParamsMPAD(Params);
  
% [G,Lam1,Var1,Len2,Var2,Mu2,vary] = unpackParamsMPAD(Params);
% Unpacks all of the parameters out of the the structure for AR2 carriers
% or
% [G,Lam1,Var1,Len2,Var2,Mu2,vary,om] = unpackParamsMPAD(Params);
% Unpacks all of the parameters out of the the structure for pFB carriers
%
% INPUTS
% Params = structure including
%  G,Lam1,Len2,Var1,Var2,Mu2,vary and optionally 'om'
%
% OUTPUTS
% G,Lam1,Len2,Var1,Var2,Mu2,vary and optionally 'om'
%
% see test_packParamsMPAD.m for tests

G = Params.G;
Lam1 = Params.Lam1;
Var1 = Params.Var1;
Len2 = Params.Len2;
Var2 = Params.Var2;
Mu2 = Params.Mu2;
vary = Params.vary;

if isfield(Params,'om')
  om = Params.om;
  varargout{1} = om; 
end
    
    