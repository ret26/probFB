function [theta,nlml] = vfeTrain(thetaInit,covfunc,x,y,noEvals,varargin)

optU = 1;
if ~isempty(varargin)
    Xu = varargin{1};
    optU = 0;
end
fname = 'vfeObjFunc';
if optU
%     d = checkgrad(fname,thetaInit,1e-7,covfunc,x,y)
%     keyboard
    [theta,f] = minimize(thetaInit,fname,-noEvals,covfunc,x,y);
else
%     d = checkgrad(fname,thetaInit,1e-7,covfunc,x,y,Xu)
%     keyboard
    [theta,f] = minimize(thetaInit,fname,-noEvals,covfunc,x,y,Xu);
end
nlml = f(end);
end
