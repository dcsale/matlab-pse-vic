clear all;
close all;
clc;

% plot styles
styles = {'k','b','r','g','m','c','--k','--b','--r','--g','--m','--c'};

%==========================================================================
%= Kernels
%==========================================================================
% polynomial coefficients (even polynomials)
p02 = [1];
p04 = [-0.5 0 2];
p06 = [1/8 0 -3/2 0 3];
p08 = [-1/48 0 1/2 0 -3 0 4];
p10 = [1/384 0 -5/48 0 5/4 0 -5 0 5];

%==========================================================================
%= Discrete moments
%==========================================================================
m = [0 2 4 6 8]; % moments
N = 128; % discretisation
sigmas = (0.1:0.1:8); % sigmas

xmin = 0;
xmax = 2;

% Allocate moment variables
moment = cell(length(m),1);

x = linspace(0,4,2*N+1); x = x(1:end-1);
dx = x(2) - x(1); dy = dx;
xc   = (xmax + xmin)/2;

[xg, yg] = meshgrid(x-xc,x-xc);
rg = sqrt((xg).^2 + (yg).^2);

%--------------------------------------------------------------------------
%- Varying sigma/dx
%--------------------------------------------------------------------------
sigmags = dx*sigmas;
for n = 1:length(sigmags)
    sigmag = sigmags(n);
    
    % Setup mesh
    rg    = sqrt((xg).^2 + (yg).^2);
    rhog  = rg/sigmag;
    rho2g = rhog.^2;
    
    % Smoothing function
    zetag{1}  = polyval(p02,rhog)/(2*pi*sigmag^2) .* exp(-rho2g/2);
    zetag{2}  = polyval(p04,rhog)/(2*pi*sigmag^2) .* exp(-rho2g/2);
    zetag{3}  = polyval(p06,rhog)/(2*pi*sigmag^2) .* exp(-rho2g/2);
    zetag{4}  = polyval(p08,rhog)/(2*pi*sigmag^2) .* exp(-rho2g/2);
    zetag{5}  = polyval(p10,rhog)/(2*pi*sigmag^2) .* exp(-rho2g/2);
    
    %----------------------------------------------------------------------
    %- Calculate the discrete moments
    %----------------------------------------------------------------------
    for i = 1:length(zetag)
        p = 1;
        for k = m
            moment{p}(i,n) = sum(sum(rg.^k.*zetag{i}))*dx*dy;
            p = p + 1;
        end
    end
end
%--------------------------------------------------------------------------
%- Plot varying sigma/dx
%--------------------------------------------------------------------------
for i = 1:length(zetag)
    figure(i+10)
    for k = 1:length(m)
        semilogy(sigmas,abs(moment{k}(i,:)),styles{k})
        hold on
    end
    title(['\zeta' num2str(2*i) ' kernel - N = ' num2str(N)])
    xlabel('q')
    ylabel('int(\rho^\beta \zeta_m)')
    %    legend('0th moment','2nd moment','4th moment','6th moment','8th moment','10th moment')
end

