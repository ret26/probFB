function [lamx,varx,dft] = freq2AR2Old(fmax,df,varMa,varargin)

% function [lamx,varx,dft] = freq2AR2Old(fmax,df,varMa,varargin)
%
% Performs a linesearch to find the parameters of an
% AR(2) process which has a specrum with centre
% frequency fmax and bandwidth df, as well as marginal
% variance varMa.
%
% INPUTS
% fmax = centre frequencies, size [D,1]
% df = bandwidths, size [D,1]
% varMa = marginal variances, size [D,1]
% OPTIONAL: verbose output (= 1 => display information)
%
% OUTPUTS
% lamx = AR(2) dynamical parameter [D,2]
% varx = AR(2) innovations noise, [D,1]
% dft = true bandwidth of the fit, [D,1]

D = length(fmax);

% Values of lam2 to check over
lam2 = linspace(-1,0,10000);

lamx = zeros(D,2);
varx = zeros(D,1);
dft = zeros(D,1);

for d=1:D

  cosom = cos(2*pi*fmax(d));
  
  % Values of lam1 that correspond to the centre frequency
  lam1 = 4*lam2./(lam2-1)*cosom;

  c1 = -8-16*cosom^2;
  c2 = -4./lam2.*(1+16*lam2.^2*cosom^2./(lam2-1).^2+lam2.^2);

  cosfd1 = cosom+1/4*sqrt(c1+c2);
  cosfd2 = cosom-1/4*sqrt(c1+c2);

  % corresponding bandwidths
  dfs = abs(acos(cosfd1)-acos(cosfd2))/(2*pi);


  % closest bandwidth to the desired bandwidth
  [val,loc] = min(abs(dfs-df(d)));

  %subplot(1,2,1)
  %hold on
  %plot(lam2,dfs)
  %plot(lam2(loc),dfs(loc),'.r','markersize',9)

  % organising output
  dft(d) = dfs(loc);
  l2 = lam2(loc);
  l1 = 4*l2/(l2-1)*cosom;
  varx(d) = varMa(d)*(1-l2-l1^2-l2^2+l2^3-l1^2*l2)/(1-l2);
  lamx(d,:) = [l1;l2]';

end

  if isempty(varargin)
    verbose=0;
  else
    verbose = varargin{1};
  end

  if verbose ==1
    for d=1:D
    str = ['Target bandwidth ', num2str(df(d)),' True bandwidth ',num2str(dft(d))];
    disp(str);
    end
  end
