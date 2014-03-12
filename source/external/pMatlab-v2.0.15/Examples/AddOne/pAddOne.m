%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple script that performs Y = X + 1.
% To run in serial without distributed arrays, set
%   PARALLEL = 0
% At the Matlab prompt type
%   pAddOne
% To run in serial with distributed arrays, set
%   PARALLEL = 1
% At the Matlab prompt type
%   pAddOne
% To run in parallel with distributed arrays
% at the Matlab prompt type
%   eval(pRUN('pAddOne',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;                  % Set matrix size.
PARALLEL = 1;    % Set control flag.
XYmap = 1;         % Create serial map.
if PARALLEL
  XYmap = map([Np 1],{},0:Np-1);   % Create parallel map.
end
X = zeros(N,N,XYmap);     % Create distributed X.
Y = zeros(N,N,XYmap);     % Create distributed Y.
Xloc = local(X);          % Get local part of X.
Yloc = local(Y);          % Get local part of Y.
Yloc = Xloc + 1;          % Add one to local part.
Y = put_local(Y,Yloc);       % Put back into Y.
