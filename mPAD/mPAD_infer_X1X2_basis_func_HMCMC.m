function [XS,acc] = mPAD_infer_X1X2_basis_func_HMCMC(x,y,varp,mup,logvars,as,bs,options,epsilon,Tau,no_samps)

% XS = HMCMC(x,variables,epsion,Tau,no_samps)
% carriers out Hamiltonian Monte Carlo for HBAD according 
% to Mackay 2003 pg 388
%
% INPUTS
% x = variables over which to sample
% variables = additional variables to be passed to the gradient and
% finction evaluation routines [gradE(x,variables) and findE(x,variables) ] 
% epsilon = size of the leap-frog step
% Tau = number of leap-frog steps
% no_samps = number of samples to return
%
% OUTPUTS
% XS = samples from the distribution
% acc = number of samples accepted

T = length(x);
g = gradE(x,y,varp,mup,logvars,as,bs,options); % set gradient using initial XS
E = findE(x,y,varp,mup,logvars,as,bs,options); % set the objective function too
samp = 1; XS = zeros(T,samp); acc = 0;

for samp=1:no_samps   % loop over samples
    p = randn(T,1);   % initialise momentum
    H = p'*p/2 + E;     % evaluate the Hamiltonian
    
    xnew = x; gnew = g;
    
    for tau = 1:Tau     % make tau leap-frog steps
        p = p - epsilon*gnew/2;         % make a half step in p
        xnew = xnew + epsilon*p;        % make a step in x
        gnew = gradE(xnew,y,varp,mup,logvars,as,bs,options);    % find new gradient
        p = p - epsilon*gnew/2;         % make a half step in p
    end
    
    Enew = findE(xnew,y,varp,mup,logvars,as,bs,options);
    Hnew = p'*p/2+Enew; % find the new value of H
    dH = Hnew - H;      % decide whether to accept
    
    if dH < 0 
        accept = 1;
    elseif rand< exp(-dH)
        accept = 1;
    else
        accept = 0;
    end
    
    if accept==1
        g = gnew; x = xnew; E = Enew;
        acc = acc+1;
    end
    XS(:,samp) = x;
    str1 = ['%% complete: ',num2str(100*samp/no_samps)];
    str2 = [' fraction accepted: ',num2str(acc/samp),'\r'];
    fprintf(str1);fprintf(str2);
end

    
        
        