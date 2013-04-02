function [X2,acc] = mPAD_infer_X1X2_basis_func_HMCMC(y,Z2,Params,opts)

% 
% carriers out Hamiltonian Monte Carlo for HBAD according 
% to Mackay 2003 pg 388
%
% INPUTS
% y = Obsevations [T,1] 
% Z2 = initial estimates for transformed envelope coefficients [T,K] 
% Params = structure containing the following parameters
%   G = modulator weights [D,K]
%   Lam1 = Carrier dynamics [D,2]
%   Var1 = Carrier noise [1,D]
%   Len2 = Modulator Length Scales [K,1] (not required)
%   Var2 = Modulator Noise [1,K]
%   vary = observation noise  
% opts = structure of options including:
%    epsilon = size of the leap-frog step
%    no_leap = number of leap-frog steps
%    no_samps = number of samples to return
%
% OUTPUTS
% X2 = one sample of transformed modulator from the distribution [T,K]
% acc = number of samples accepted

tol = 7;

[G,Lam1,Var1,Len2,Var2,Mu2,vary,om] = unpackParamsMPAD(Params);

[T,K] = size(Z2);
z2 = Z2(:);

[Obj,dObj] = getObj_mPAD_noG_basis_func(z2,y,Params,tol);
g = dObj*T; % set gradient using initial z2
E = Obj*T; % set the objective function too

acc = 0;

figure
for samp=1:opts.no_samps   % loop over samples
    p = randn(T*K,1);   % initialise momentum
    H = p'*p/2 + E;     % evaluate the Hamiltonian
    
    z2new = z2; gnew = g;
    
    for tau = 1:opts.no_leap     % make no_leap leap-frog steps
        p = p - opts.epsilon*gnew/2;         % make a half step in p
        z2new = z2new + opts.epsilon*p;        % make a step in z2
        
	[ObjNew,dObjNew] = getObj_mPAD_noG_basis_func(z2new,y,Params,tol);
	gnew = dObjNew*T; % find new gradient
	Enew = ObjNew*T;
        p = p - opts.epsilon*gnew/2;         % make a half step in p
    end
        
    Hnew = p'*p/2+Enew; % find the new value of H
    dH = Hnew - H;      % decide whether to accept
    dH
    if dH < 0 
        accept = 1;
    elseif rand< exp(-dH)
        accept = 1;
    else
        accept = 0;
    end
    
    if accept==1
        g = gnew; z2 = z2new; E = Enew;
        acc = acc+1;
    end
%    XS(:,samp) = x;
    str1 = ['%% complete: ',num2str(100*samp/opts.no_samps)];
    str2 = [' fraction accepted: ',num2str(acc/samp),'\r'];
    fprintf(str1);fprintf(str2);

    Z2 = reshape(z2,[T,K]);      
    X2 = Z2_to_X2(Z2,Len2,Var2,tol);

    plot(X2)
    drawnow
    
end

Z2 = reshape(z2,[T,K]);      
X2 = Z2_to_X2(Z2,Len2,Var2,tol);
       
        