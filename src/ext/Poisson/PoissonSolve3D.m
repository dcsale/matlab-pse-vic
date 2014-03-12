clear all; close all; clc;
%**************************************************************************
%*
%* Program:      PoissonSolve3D.m
%* 
%* Description:  A 3D high order converging Poisson solver for unbounded 
%*               domains based on a Greens function solution. The method
%*               is applied to the bump function test case.
%*                
%* Publication:  'A high order solver for the unbounded Poisson equation'
%*               J. Comput. Phys.
%*
%* Authors:      Mads Mølholm Hejlesen (a)
%*               Johannes Tophøj Rasmussen (a)
%*               Philippe Chatelain (b)
%*               Jens Honore Walther (a,c)*
%*
%*               (a) Technical University of Denmark
%*               (b) Universite catholique de Louvain
%*               (c) ETH Zürich
%*               *Corresponding author at jhw@mek.dtu.dk
%*
%**************************************************************************

%==========================================================================
% Input options 
%==========================================================================
% Regularisation order of integration kernel: 2,4,...,10 (0 = Non-regularised)
kernel    = 4;

% Solution
% 0 = solve for stream function by G kernel
% 1 = solve for velocity by K kernels, 
% 2 = solve for velocity by G kernel + spectral differentiating
solve_vel = 2; 

% Smoothing radius relative to mesh size: epsilon = alpha*dx (default 2)
alpha     = 2; 

% Number of mesh cells (use vector for convergence studies)
NXs       = 32*2.^(0:2);  

% Testcases
% 1 = Bump function: SPHERICAL SCALAR FIELD (only for solve_vel = 0)
% 2 = Bump function: VORTEX RING (xy-plane)
testcase = 2;

%==========================================================================
%==========================================================================
%--------------------------------------------------------------------------
%- Setup domain range and length
%--------------------------------------------------------------------------
if(testcase == 1)
    xmin = [-1 -1 -1]; xmax = [1 1 1];
elseif(testcase == 2)
    xmin = [0 0 0]; xmax = [2 2 1];
end
L = xmax-xmin; % Domain length

%--------------------------------------------------------------------------
%- Loop over discretisation levels
%--------------------------------------------------------------------------
for N = 1:length(NXs)
    % Number of grid points in each direction
    NX(1) = NXs(N);
    NX(2) = round(L(2)/L(1)*NX(1));
    NX(3) = round(L(3)/L(1)*NX(1));
    %----------------------------------------------------------------------
    % Setup domain (_ext is the extended domain)
    %----------------------------------------------------------------------
    x{1} = linspace(xmin(1),xmax(1),NX(1));
    x{2} = linspace(xmin(2),xmax(2),NX(2));
    x{3} = linspace(xmin(3),xmax(3),NX(3));
    
    [fieldx, fieldy, fieldz] = ndgrid(x{1},x{2},x{3});
    
    % Grid spacing
    dx(1) = x{1}(2)-x{1}(1);
    dx(2) = x{2}(2)-x{2}(1);
    dx(3) = x{3}(2)-x{3}(1);
    
    % Centred domain
    centre = (xmin+xmax)/2;
    fieldx_cen{1} = fieldx - centre(1);
    fieldx_cen{2} = fieldy - centre(2);
    fieldx_cen{3} = fieldz - centre(3);
    
    % Extended domain (FFT shifted)
    x_ext{1} = x{1}-xmin(1); x_ext{1} = [x_ext{1} -x_ext{1}(end)-dx(1) -x_ext{1}(end:-1:2)];
    x_ext{2} = x{2}-xmin(2); x_ext{2} = [x_ext{2} -x_ext{2}(end)-dx(2) -x_ext{2}(end:-1:2)];
    x_ext{3} = x{3}-xmin(3); x_ext{3} = [x_ext{3} -x_ext{3}(end)-dx(3) -x_ext{3}(end:-1:2)];
    
    [fieldx_ext{1}, fieldx_ext{2}, fieldx_ext{3}] = ndgrid(x_ext{1},x_ext{2},x_ext{3});
    fieldr_ext = sqrt(fieldx_ext{1}.^2 + fieldx_ext{2}.^2 + fieldx_ext{3}.^2);
    
    if(solve_vel == 2)
        % Wavenumbers for spectral differentiating
        ks = 1/dx(1);
        k_ext{1} = 2*pi*linspace(-ks/2,ks/2,2*NX(1)+1);
        k_ext{1} = fftshift(k_ext{1}(1:end-1));
        
        ks = 1/dx(2);
        k_ext{2} = 2*pi*linspace(-ks/2,ks/2,2*NX(2)+1);
        k_ext{2} = fftshift(k_ext{2}(1:end-1));
        
        ks = 1/dx(3);
        k_ext{3} = 2*pi*linspace(-ks/2,ks/2,2*NX(3)+1);
        k_ext{3} = fftshift(k_ext{3}(1:end-1));
        
        [fieldk{1}, fieldk{2}, fieldk{3}] = ndgrid(k_ext{1},k_ext{2},k_ext{3});
    end
    
    clear fieldx fieldy fieldz
    %----------------------------------------------------------------------
    % Vorticity distribution
    %----------------------------------------------------------------------
    vort{1} = zeros(NX(1),NX(2),NX(3));
    vort{2} = zeros(NX(1),NX(2),NX(3));
    vort{3} = zeros(NX(1),NX(2),NX(3));
    
    if(testcase == 1)
        % Bump function: SPHERICAL SCALAR FIELD
        c = 20; R = 1; % Function constants
        fieldr_cen = sqrt(fieldx_cen{1}.^2 + fieldx_cen{2}.^2 + fieldx_cen{3}.^2);
        
        vort{2}(fieldr_cen < R) = 2*c*R^2*exp(-c*R^2./((R - fieldr_cen(fieldr_cen < R))...
            .*(R + fieldr_cen(fieldr_cen < R)))).*(3*R^4 - 2*R^2*fieldr_cen(fieldr_cen < R).^2 ...
            - fieldr_cen(fieldr_cen < R).^4 - 2*c*R^2*fieldr_cen(fieldr_cen < R).^2) ...
            .*(R^2 - fieldr_cen(fieldr_cen < R).^2).^(-4);
    elseif(testcase == 2)
        % Bump function: VORTEX RING (in the xy-plane)
        c = 10; R = 0.5; % Function constants
        vort_mag      = zeros(NX(1),NX(2),NX(3));
        fieldr_cen    = sqrt(fieldx_cen{1}.^2 + fieldx_cen{2}.^2);
        fieldphi_cen  = sqrt((fieldr_cen-R).^2 + fieldx_cen{3}.^2);
        fieldtheta    = atan2(fieldx_cen{2},fieldx_cen{1});
        
        vort_mag(fieldphi_cen < R) = exp(-c*R^2./(2*R*fieldr_cen(fieldphi_cen < R) - ...
            fieldr_cen(fieldphi_cen < R).^2 - fieldx_cen{3}(fieldphi_cen < R).^2)) ...
            .*(4*c^2*R^4*fieldx_cen{3}(fieldphi_cen < R).^2.*fieldr_cen(fieldphi_cen < R).^2 ...
            - 16*R^4*fieldr_cen(fieldphi_cen < R).^4 ...
            + 32*R^3*fieldr_cen(fieldphi_cen < R).^5 ...
            - 24*R^2*fieldr_cen(fieldphi_cen < R).^6 ...
            + 8*R*fieldr_cen(fieldphi_cen < R).^7 ...
            - 4*fieldr_cen(fieldphi_cen < R).^6.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
            - 6*fieldr_cen(fieldphi_cen < R).^4.*fieldx_cen{3}(fieldphi_cen < R).^4 ...
            - 4*fieldr_cen(fieldphi_cen < R).^2.*fieldx_cen{3}(fieldphi_cen < R).^6 ...
            - 8*c*R^5*fieldr_cen(fieldphi_cen < R).^3 ...
            + 8*c*R^4*fieldr_cen(fieldphi_cen < R).^4 ...
            - 6*c*R^3*fieldr_cen(fieldphi_cen < R).^5 ...
            + 4*c^2*R^6*fieldr_cen(fieldphi_cen < R).^2 ...
            - 8*c^2*R^5*fieldr_cen(fieldphi_cen < R).^3 ...
            + 4*c^2*R^4*fieldr_cen(fieldphi_cen < R).^4 ...
            + 2*c*R^2*fieldr_cen(fieldphi_cen < R).^6 ...
            + 32*R^3*fieldr_cen(fieldphi_cen < R).^3.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
            - 48*R^2*fieldr_cen(fieldphi_cen < R).^4.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
            - 24*R^2*fieldr_cen(fieldphi_cen < R).^2.*fieldx_cen{3}(fieldphi_cen < R).^4 ...
            + 24*R*fieldr_cen(fieldphi_cen < R).^5.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
            + 24*R*fieldr_cen(fieldphi_cen < R).^3.*fieldx_cen{3}(fieldphi_cen < R).^4 ...
            + 8*R*fieldr_cen(fieldphi_cen < R).*fieldx_cen{3}(fieldphi_cen < R).^6 ...
            + 2*c*R^3*fieldr_cen(fieldphi_cen < R).*fieldx_cen{3}(fieldphi_cen < R).^4 ...
            + 2*c*R^2*fieldr_cen(fieldphi_cen < R).^2.*fieldx_cen{3}(fieldphi_cen < R).^4 ...
            - 4*c*R^3*fieldr_cen(fieldphi_cen < R).^3.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
            + 4*c*R^2*fieldr_cen(fieldphi_cen < R).^4.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
            - fieldr_cen(fieldphi_cen < R).^8 - fieldx_cen{3}(fieldphi_cen < R).^8) ...
            .*(2*R*fieldr_cen(fieldphi_cen < R) - ...
            fieldr_cen(fieldphi_cen < R).^2 - fieldx_cen{3}(fieldphi_cen < R).^2).^(-4) ...
            .*fieldr_cen(fieldphi_cen < R).^(-2);
        
        vort{1} = sin(fieldtheta).*vort_mag;
        vort{2} = -cos(fieldtheta).*vort_mag;
        vort{3} = zeros(NX(1),NX(2),NX(3));
    end
    
    %----------------------------------------------------------------------
    % Analytical solution
    %----------------------------------------------------------------------
    if(solve_vel == 0)
        disp(['Solving stream function using kernel G' num2str(kernel) ' at NX = ' num2str(NX)])
        stream_ana{1} = zeros(NX(1),NX(2),NX(3));
        stream_ana{2} = zeros(NX(1),NX(2),NX(3));
        stream_ana{3} = zeros(NX(1),NX(2),NX(3));
        
        if(testcase == 1)
            % Bump function: SPHERICAL SCALAR FIELD
            stream_ana{2}(fieldr_cen < R) = exp(-c./(1 - fieldr_cen(fieldr_cen < R).^2/R^2));
            
        elseif(testcase == 2)
            % Bump function: VORTEX RING (in the xy-plane)
            stream_ana_mag = zeros(NX(1),NX(2),NX(3));
            stream_ana_mag(fieldphi_cen < R) = exp(-c./(1 - (fieldr_cen(fieldphi_cen < R).^2 ...
                + R^2 + fieldx_cen{3}(fieldphi_cen < R).^2 - 2*R*fieldr_cen(fieldphi_cen < R))/R^2));
            stream_ana{1} = -sin(fieldtheta).*stream_ana_mag;
            stream_ana{2} = cos(fieldtheta).*stream_ana_mag;
            stream_ana{3} = zeros(NX(1),NX(2),NX(3));
        end
    else
        if(solve_vel == 1)
            disp(['Solving velocity using kernel K' num2str(kernel) ' at NX = ' num2str(NX)])
        else
            disp(['Solving velocity using kernel G' num2str(kernel) ' at NX = ' num2str(NX)])
        end
        vel_ana{1} = zeros(NX(1),NX(2),NX(3));
        vel_ana{2} = zeros(NX(1),NX(2),NX(3));
        vel_ana{3} = zeros(NX(1),NX(2),NX(3));
        
        if(testcase == 1)
            % Bump function: SPHERICAL SCALAR FIELD
            error('Velocity field is not available for this test function - set solve_vel = 0')
        elseif(testcase == 2)
            % Bump function VORTEX RING (in the xy-plane)
            vel_ana_r = zeros(NX(1),NX(2),NX(3));
            vel_ana_r(fieldphi_cen < R) = 2*c*R^2*fieldx_cen{3}(fieldphi_cen < R) ...
                .* exp(-c*R^2./(2*R*fieldr_cen(fieldphi_cen < R) - ...
                fieldr_cen(fieldphi_cen < R).^2 - fieldx_cen{3}(fieldphi_cen < R).^2)) ...
                .*(2*R*fieldr_cen(fieldphi_cen < R) - fieldr_cen(fieldphi_cen < R).^2 ...
                - fieldx_cen{3}(fieldphi_cen < R).^2).^(-2);
            
            vel_ana{1} = cos(fieldtheta).*vel_ana_r;
            vel_ana{2} = sin(fieldtheta).*vel_ana_r;
            
            vel_ana{3}(fieldphi_cen < R) = exp(-c*R^2./(2*R*fieldr_cen(fieldphi_cen < R) - ...
                fieldr_cen(fieldphi_cen < R).^2 - fieldx_cen{3}(fieldphi_cen < R).^2)) ...
                .*(4*R^2*fieldr_cen(fieldphi_cen < R).^2 ...
                - 4*R*fieldr_cen(fieldphi_cen < R).^3 ...
                - 4*R*fieldr_cen(fieldphi_cen < R).*fieldx_cen{3}(fieldphi_cen < R).^2 ...
                + fieldr_cen(fieldphi_cen < R).^4 ...
                + 2*fieldr_cen(fieldphi_cen < R).^2.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
                + fieldx_cen{3}(fieldphi_cen < R).^4 ...
                + 2*c*R^3*fieldr_cen(fieldphi_cen < R) ...
                - 2*c*R^2*fieldr_cen(fieldphi_cen < R).^2) ...
                .*(2*R*fieldr_cen(fieldphi_cen < R) - fieldr_cen(fieldphi_cen < R).^2 ...
                - fieldx_cen{3}(fieldphi_cen < R).^2).^(-2).*fieldr_cen(fieldphi_cen < R).^(-1);
        end        
    end
    
    clear fieldx_cen{1} fieldr_cen
    
    %######################################################################
    % THE POISSON SOLVER
    %######################################################################
    %----------------------------------------------------------------------
    % Fourier transform vorticity field
    %----------------------------------------------------------------------
    for i = 1:3
        vort_ext{i} = zeros(2*NX(1),2*NX(2),2*NX(3));
        vort_ext{i}(1:NX(1),1:NX(2),1:NX(3)) = vort{i};
        vort_fft{i} = fftn(vort_ext{i});
    end
    
    clear vort_ext
    %----------------------------------------------------------------------
    % Integration kernel
    %----------------------------------------------------------------------
    epsilon = alpha*max(dx); % smoothing radius
    rho   = fieldr_ext/epsilon; % normalised radial distance
    if(kernel == 0) % 0th order (No regularisation)
        if(solve_vel == 1)
            K{1} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{1};
            K{1}(1,1,1) = 0;
            K{2} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{2};
            K{2}(1,1,1) = 0;
            K{3} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{3};
            K{3}(1,1,1) = 0;
        else
            G = 1./(4*pi*fieldr_ext);
            G(1,1,1) = 1;
        end
    elseif(kernel == 2) % 2nd order
        if(solve_vel == 1)
            K{1} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{1}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{1}(1,1,1) = 0;
            K{2} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{2};
            K{2}(1,1,1) = 0;
            K{3} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{3};
            K{3}(1,1,1) = 0;            
        else
            G = 1./(4*pi*fieldr_ext).*erf(rho/sqrt(2));
            G(1,1,1) = sqrt(2)/(4*pi^(3/2)*epsilon);
        end
    elseif(kernel == 4) % 4th order
        if(solve_vel == 1)
            K{1} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{1}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{1}(1,1,1) = 0;
            K{2} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{2}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{2}(1,1,1) = 0;
            K{3} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{3}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{3}(1,1,1) = 0;
        else
            G = 1./(4*pi*fieldr_ext).*((rho).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
            G(1,1,1) = 3*sqrt(2)/(8*pi^(3/2)*epsilon);
        end
    elseif(kernel == 6) % 6th order
        if(solve_vel == 1)
            K{1} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{1}.*((-2*rho + 9/4*rho.^3 - 1/4*rho.^5).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{1}(1,1,1) = 0;
            K{2} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{2}.*((-2*rho + 9/4*rho.^3 - 1/4*rho.^5).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{2}(1,1,1) = 0;
            K{3} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{3}.*((-2*rho + 9/4*rho.^3 - 1/4*rho.^5).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{3}(1,1,1) = 0;            
        else
            G = 1./(4*pi*fieldr_ext).*((7/4*rho - 1/4*rho.^3).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
            G(1,1,1) = 15*sqrt(2)/(32*pi^(3/2)*epsilon);
        end
    elseif(kernel == 8) % 8th order
        if(solve_vel == 1)
            K{1} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{1}.*((-2*rho + 89/24*rho.^3 - 20/24*rho.^5 + 1/24*rho.^7).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{1}(1,1,1) = 0;
            K{2} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{2}.*((-2*rho + 89/24*rho.^3 - 20/24*rho.^5 + 1/24*rho.^7).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{2}(1,1,1) = 0;
            K{3} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{3}.*((-2*rho + 89/24*rho.^3 - 20/24*rho.^5 + 1/24*rho.^7).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{3}(1,1,1) = 0;            
        else
            G = 1./(4*pi*fieldr_ext).*((19/8*rho - 2/3*rho.^3 + 1/24*rho.^5).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
            G(1,1,1) = 35*sqrt(2)/(64*pi^(3/2)*epsilon);
        end
    elseif(kernel == 10) % 10th order
        if(solve_vel == 1)
            K{1} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{1}.*((-2*rho + 1027/192*rho.^3 - 349/192*rho.^5 + 35/192*rho.^7 - 1/192*rho.^9).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{1}(1,1,1) = 0;
            K{2} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{2}.*((-2*rho + 1027/192*rho.^3 - 349/192*rho.^5 + 35/192*rho.^7 - 1/192*rho.^9).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{2}(1,1,1) = 0;
            K{3} = - 1./(4*pi*fieldr_ext.^3).*fieldx_ext{3}.*((-2*rho + 1027/192*rho.^3 - 349/192*rho.^5 + 35/192*rho.^7 - 1/192*rho.^9).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
            K{3}(1,1,1) = 0;            
        else
            G = 1./(4*pi*fieldr_ext).*((187/64*rho - 233/192*rho.^3 + 29/192*rho.^5 - 1/192*rho.^7).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
            G(1,1,1) = 315*sqrt(2)/(512*pi^(3/2)*epsilon);
        end      
    else
        error('Specified kernel is not implemented')
    end
    
    % Fourier transform integration kernel (normalised)
    if(solve_vel == 1)
        K_fft{1} = fftn(K{1}*dx(1)*dx(2)*dx(3));
        K_fft{2} = fftn(K{2}*dx(1)*dx(2)*dx(3));
        K_fft{3} = fftn(K{3}*dx(1)*dx(2)*dx(3));
    else
        G_fft = fftn(G*dx(1)*dx(2)*dx(3));
    end

    %----------------------------------------------------------------------
    % Convolve and calculate the curl in Fourier space
    %----------------------------------------------------------------------
    if(solve_vel == 1)
        vel_fft{1} = vort_fft{3}.*K_fft{2} - vort_fft{2}.*K_fft{3};
        vel_fft{2} = vort_fft{1}.*K_fft{3} - vort_fft{3}.*K_fft{1};
        vel_fft{3} = vort_fft{2}.*K_fft{1} - vort_fft{1}.*K_fft{2};
    else
        stream_fft{1} = vort_fft{1}.*G_fft;
        stream_fft{2} = vort_fft{2}.*G_fft;
        stream_fft{3} = vort_fft{3}.*G_fft;
        if(solve_vel == 2)
            vel_fft{1} = stream_fft{3}.*(sqrt(-1)*fieldk{2}) - stream_fft{2}.*(sqrt(-1)*fieldk{3});
            vel_fft{2} = stream_fft{1}.*(sqrt(-1)*fieldk{3}) - stream_fft{3}.*(sqrt(-1)*fieldk{1});
            vel_fft{3} = stream_fft{2}.*(sqrt(-1)*fieldk{1}) - stream_fft{1}.*(sqrt(-1)*fieldk{2});
        end
    end
    
    %----------------------------------------------------------------------
    % Transform back to physical space
    %----------------------------------------------------------------------
    if(solve_vel == 0)
        stream_ext{1} = real(ifftn(stream_fft{1}));
        stream_ext{2} = real(ifftn(stream_fft{2}));
        stream_ext{3} = real(ifftn(stream_fft{3}));
    else
        vel_ext{1} = real(ifftn(vel_fft{1}));
        vel_ext{2} = real(ifftn(vel_fft{2}));
        vel_ext{3} = real(ifftn(vel_fft{3}));
    end
    
    %----------------------------------------------------------------------
    % Extract numerical solution
    %----------------------------------------------------------------------
    if(solve_vel == 0)
        stream{1} = stream_ext{1}(1:NX(1),1:NX(2),1:NX(3));
        stream{2} = stream_ext{2}(1:NX(1),1:NX(2),1:NX(3));
        stream{3} = stream_ext{3}(1:NX(1),1:NX(2),1:NX(3));
        
        % Calculate diffference from analytical solution
        diff = sqrt((stream{1} - stream_ana{1}).^2 + (stream{2} ...
            - stream_ana{2}).^2 + (stream{3} - stream_ana{3}).^2) ...
            /max(max(max(sqrt(stream_ana{1}.^2 + stream_ana{2}.^2 + stream_ana{3}.^2))));
    else
        vel{1} = vel_ext{1}(1:NX(1),1:NX(2),1:NX(3));
        vel{2} = vel_ext{2}(1:NX(1),1:NX(2),1:NX(3));
        vel{3} = vel_ext{3}(1:NX(1),1:NX(2),1:NX(3));
        
        % Calculate diffference from analytical solution
        diff = sqrt((vel{1} - vel_ana{1}).^2 + (vel{2} ...
            - vel_ana{2}).^2 + (vel{3} - vel_ana{3}).^2) ...
            /max(max(max(sqrt(vel_ana{1}.^2 + vel_ana{2}.^2 + vel_ana{3}.^2))));
    end
    
    %----------------------------------------------------------------------
    % Error calculation (rms, max and L2-norm)
    %----------------------------------------------------------------------
    err_rms(N) = sqrt(sum(sum(sum(diff.^2)))/(NX(1)*NX(2)*NX(3)));
    err_max(N) = max(max(max(abs(diff))));    
    err_L2(N)  = sum(sum(sum(diff.^2)))^(1/2)/(NX(1)*NX(2)*NX(3));
    
    dxs(N)     = dx(1);
end
%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
if(length(NXs) > 1)
    % Convergence test
    figure(10)
    set(gca,'FontSize',14)
    loglog(dxs,err_L2,'b')
    hold on
    loglog(dxs,err_max,'g')
    loglog(dxs,err_rms,'r')
    loglog(dxs,dxs.^2,'k')
    loglog(dxs,dxs.^4,'k')
    loglog(dxs,dxs.^6,'k')
    loglog(dxs,dxs.^8,'k')
    loglog(dxs,dxs.^10,'k')
    xlim([min(dxs)/2 max(dxs)*2])
    title('Convergence test')
    xlabel('dx')
    ylabel('error')
    legend('rms error','max error','L2 error','dx^2','dx^4','dx^6','dx^8','dx^{10}',1)
end

% Initial and resulting field along (x,y,z) = (:,0,0)
figure(20)
set(gca,'FontSize',14)
plot(x{1},vort{2}(:,NX(2)/2,NX(3)/2)/max(max(max(vort{2}))),'b')
hold on
if(solve_vel == 0)
    plot(x{1},stream_ana{2}(:,NX(2)/2,NX(3)/2)/max(max(max(stream_ana{2}))),'r')
    plot(x{1},stream{2}(:,NX(2)/2,NX(3)/2)/max(max(max(stream{2}))),'--g')
    title('Initial vorticity and resulting stream field at (x,y,z) = (:,0,0)')
    xlabel('x')
    ylabel('normalised \omega_2 , \psi_3')
    legend('vorticity','analytic stream function','calculated stream function',4)
else
    plot(x{1},vel_ana{3}(:,NX(2)/2,NX(3)/2)/max(max(max(vel_ana{3}))),'r')
    plot(x{1},vel{3}(:,NX(2)/2,NX(3)/2)/max(max(max(vel{3}))),'--g')  
    title('Initial vorticity and resulting velocity field at (x,y,z) = (:,0,0)')
    xlabel('x')
    ylabel('normalised \omega_2 , u_3')
    legend('vorticity','analytic velocity','calculated velocity',4)
end

