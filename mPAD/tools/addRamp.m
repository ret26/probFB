function Yramp = addRamp(Y,lens)

% function Yramp = addRamp(Y,lens)
%
% Adds a cosine ramp to the start and the end of a stimulus. Does
% this down the columns of Y.
%
% INPUTS
% Y = data [T,D]
% lens = length-scale of cosine ramp [D,1]
%
% OUTPUTS
% Yramp = ramps added to the end and start of Y [T,D]

[T,D] = size(Y);
lens = ceil(lens);
Yramp = Y;

for d=1:D
  onset = sin(2*pi*[1:lens(d)]'/(4*lens(d)));
  Yramp(1:lens(d),d) = Yramp(1:lens(d),d).*onset;

  Yramp(T-lens(d)+1:T,d) = Yramp(T-lens(d)+1:T,d).*onset(end:-1:1);
end
