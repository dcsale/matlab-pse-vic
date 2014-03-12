clear all; close all; clc;
%**************************************************************************
%*
%* Program:      PoissonSolve2D.m
%*
%* Description:  A 2D high order converging Poisson solver for unbounded
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
% Regularisation order of integration kernel: 2,4,...,10, 100 (0 = Non-regularised)
kernel    = 100;

colour = 'b';

% Solution
% 0 = solve for stream function by G kernel
% 1 = solve for velocity by K kernels,
% 2 = solve for velocity by G kernel + spectral differentiating
solve_vel = 1;

% Smoothing radius relative to mesh size: epsilon = alpha*dx (default 2)
alpha     = 2;

% Number of grid points (use vector for convergence studies)
NXs       = 16*2.^(0:4);

% Testcases
% 1 = Solid body rotation (velocity only)
% 2 = Beale patch [2] (velocity only)
% 3 = Perlman vorticity patch [1] (velocity only)
% 4 = Bump function
testcase = 4;

%==========================================================================
%==========================================================================
%--------------------------------------------------------------------------
%- Setup domain range and length
%--------------------------------------------------------------------------
if(testcase == 1)
    xmin = [-2 -2]; xmax = [2 2];
elseif(testcase == 2)
    xmin = [-2 -2]; xmax = [2 2];
elseif(testcase == 3)
    xmin = [-1 -1]; xmax = [1 1];
elseif(testcase == 4)
    xmin = [-1 -1]; xmax = [1 1];       
end
L = xmax-xmin; % Domain length

%--------------------------------------------------------------------------
%- Loop over discretisation levels
%--------------------------------------------------------------------------
for N = 1:length(NXs)
    NX(1) = NXs(N);
    NX(2) = round(L(2)/L(1)*NX(1));
    %----------------------------------------------------------------------
    % Setup domain (_ext is the extended domain)
    %----------------------------------------------------------------------
    x{1} = linspace(xmin(1),xmax(1),NX(1));
    x{2} = linspace(xmin(2),xmax(2),NX(2));
    
    [fieldx{1} fieldx{2}] = ndgrid(x{1},x{2});
    
    % Grid spacing
    dx(1) = x{1}(2)-x{1}(1);
    dx(2) = x{2}(2)-x{2}(1);
    
    % Centred domain
    centre = (xmin+xmax)/2;
    fieldx_cen{1} = fieldx{1} - centre(1);
    fieldx_cen{2} = fieldx{2} - centre(2);
    fieldr_cen = sqrt(fieldx_cen{1}.^2 + fieldx_cen{2}.^2);
    
    % Extended domain (FFT shifted)
    x_ext{1} = x{1} - xmin(1); x_ext{1} = [x_ext{1} -x_ext{1}(end)-dx(1) -x_ext{1}(end:-1:2)];
    x_ext{2} = x{2} - xmin(2); x_ext{2} = [x_ext{2} -x_ext{2}(end)-dx(2) -x_ext{2}(end:-1:2)];
    
    [fieldx_ext{1} fieldx_ext{2}] = ndgrid(x_ext{1},x_ext{1});
    fieldr_ext = sqrt(fieldx_ext{1}.^2 + fieldx_ext{2}.^2);
    
    if(solve_vel == 2)
        % Wavenumbers for spectral differentiating
        ks = 1/dx(1);
        k_ext{1} = 2*pi*linspace(-ks/2,ks/2,2*NX(1)+1);
        k_ext{1} = fftshift(k_ext{1}(1:end-1));
        
        ks = 1/dx(2);
        k_ext{2} = 2*pi*linspace(-ks/2,ks/2,2*NX(2)+1);
        k_ext{2} = fftshift(k_ext{2}(1:end-1));
        
        [fieldk{1} fieldk{2}] = ndgrid(k_ext{1},k_ext{2});
    end
    
    %----------------------------------------------------------------------
    % Vorticity distribution
    %----------------------------------------------------------------------
    vort = zeros(NX(1),NX(2));
    if(testcase == 1)
        % Solid body rotation
        vort(fieldr_cen <= 1) = 1;
        vort(fieldr_cen > 1)  = 0;  
     elseif(testcase == 2)       
        % Beale vorticity patch
        vort(fieldr_cen <= 1) = (1-fieldr_cen(fieldr_cen <= 1).^2).^3;
        vort(fieldr_cen > 1)  = 0;        
     elseif(testcase == 3)       
        % Perlman vorticity patch
        vort(fieldr_cen <= 1) = (1-fieldr_cen(fieldr_cen <= 1).^2).^7;
        vort(fieldr_cen > 1)  = 0;
    elseif(testcase == 4)
        % Bump function
        c = 20;
        vort(fieldr_cen < 1) = -4*c*exp(-c./(1-fieldx_cen{1}(fieldr_cen < 1).^2 - fieldx_cen{2}(fieldr_cen < 1).^2)).*...
            (c*fieldx_cen{1}(fieldr_cen < 1).^2 + fieldx_cen{1}(fieldr_cen < 1).^4 + 2*fieldx_cen{1}(fieldr_cen < 1).^2.*fieldx_cen{2}(fieldr_cen < 1).^2 + c*fieldx_cen{2}(fieldr_cen < 1).^2 + fieldx_cen{2}(fieldr_cen < 1).^4 - 1).* ...
            (-1 + fieldx_cen{1}(fieldr_cen < 1).^2 + fieldx_cen{2}(fieldr_cen < 1).^2).^(-4);
        vort(fieldr_cen >= 1)  = 0;
    else        
        error('testcase unknown')
    end
    %----------------------------------------------------------------------
    % Analytic solution
    %----------------------------------------------------------------------
    if(solve_vel == 0)
        disp(['Solving stream function using kernel G' num2str(kernel) ' at NX = ' num2str(NXs(N))])
        stream_ana = zeros(NX(1),NX(2));
        if(testcase == 1)
            % Solid body patch
            error('There exists no stream function for the Solid body testcase')
        elseif(testcase == 2)            
            % Beale vorticity patch
            error('There exists no stream function for the Beale testcase')            
        elseif(testcase == 3)            
            % Perlman vorticity patch
            error('There exists no stream function for the Perlman testcase')
        elseif(testcase == 4)
            % Bump function
            stream_ana(fieldr_cen < 1) = exp(-c./(1 - fieldx_cen{1}(fieldr_cen < 1).^2 - fieldx_cen{2}(fieldr_cen < 1).^2));
        end
    else
        if(solve_vel == 1)
            disp(['Solving velocity using kernel K' num2str(kernel) ' at NX = ' num2str(NXs(N))])
        else
            disp(['Solving velocity using kernel G' num2str(kernel) ' at NX = ' num2str(NXs(N))])
        end
        vel_ana{1} = zeros(NX(1),NX(2));
        vel_ana{2} = zeros(NX(1),NX(2));
        if(testcase == 1)
            % Solid body rotation
            vel_ana{1} = -1./(2*fieldr_cen.^2).*fieldx_cen{2}.* (fieldr_cen > 1) ...
                - 0.5*fieldx_cen{2}.* (fieldr_cen <= 1);
            vel_ana{2} =  1./(2*fieldr_cen.^2).*fieldx_cen{1}.* (fieldr_cen > 1)...
                + 0.5*fieldx_cen{1}.* (fieldr_cen <= 1);
        elseif(testcase == 2)     
            % Beale vorticity patch            
            vel_ana{1} = (fieldr_cen.^6 - 4*fieldr_cen.^4 + 6*fieldr_cen.^2 - 4).*fieldx_cen{2}*(1/8).* (fieldr_cen <= 1) ...
                - 1./(8*fieldr_cen.^2).*fieldx_cen{2} .* (fieldr_cen > 1); 
            vel_ana{2} = -(fieldr_cen.^6 - 4*fieldr_cen.^4 + 6*fieldr_cen.^2 - 4).*fieldx_cen{1}*(1/8).* (fieldr_cen <= 1) ...
                + 1./(8*fieldr_cen.^2).*fieldx_cen{1} .* (fieldr_cen > 1);  
            vel_ana{1}(fieldr_cen == 0) = 0;
            vel_ana{2}(fieldr_cen == 0) = 0;
        elseif(testcase == 3)     
            % Perlman vorticity patch
            vel_ana{1} = -(1-(1-fieldr_cen.^2).^8)./(16*fieldr_cen.^2).* fieldx_cen{2} .* (fieldr_cen <= 1) ...
                - 1./(16*fieldr_cen.^2).*fieldx_cen{2} .* (fieldr_cen > 1);
            vel_ana{2} = (1-(1-fieldr_cen.^2).^8)./(16*fieldr_cen.^2).*fieldx_cen{1} .* (fieldr_cen <= 1) ...
                + 1./(16*fieldr_cen.^2).*fieldx_cen{1} .* (fieldr_cen > 1);
            vel_ana{1}(fieldr_cen == 0) = 0;
            vel_ana{2}(fieldr_cen == 0) = 0;
        elseif(testcase == 4)
            % Bump function
            vel_ana{1}(fieldr_cen < 1) = -2*c*fieldx_cen{2}(fieldr_cen < 1).*exp(-c./(1 - fieldx_cen{1}(fieldr_cen < 1).^2 - fieldx_cen{2}(fieldr_cen < 1).^2)) ...
                ./(1 - fieldx_cen{1}(fieldr_cen < 1).^2 - fieldx_cen{2}(fieldr_cen < 1).^2).^2;
            vel_ana{2}(fieldr_cen < 1) = 2*c*fieldx_cen{1}(fieldr_cen < 1).*exp(-c./(1 - fieldx_cen{1}(fieldr_cen < 1).^2 - fieldx_cen{2}(fieldr_cen < 1).^2)) ...
                ./(1 - fieldx_cen{1}(fieldr_cen < 1).^2 - fieldx_cen{2}(fieldr_cen < 1).^2).^2;
            vel_ana{1}(fieldr_cen == 0) = 0;
            vel_ana{2}(fieldr_cen == 0) = 0;
        end        
    end
    
    clear fieldx fieldx_cen fieldr_cen
    
    %##########################################################################
    % THE SOLVER
    %##########################################################################
    %----------------------------------------------------------------------
    % Fourier transform vorticity field
    %----------------------------------------------------------------------
    vort_ext = zeros(2*NX(1),2*NX(2));
    vort_ext(1:NX(1),1:NX(2)) = vort;
    vort_fft = fft2(vort_ext);
    
    clear vort_ext
    %----------------------------------------------------------------------
    % Integration kernel
    %----------------------------------------------------------------------
    epsilon = alpha*max(dx);
    gamma = 0.5772156649;
    if(kernel == 0) % 0th order (No regularisation)
        if(solve_vel == 1)
            K{1} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{1};
            K{1}(1,1) = 0;
            K{2} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{2};
            K{2}(1,1) = 0;
        else
            G = -1/(2*pi)*log(fieldr_ext);
            G(1,1) = 1;
        end
    elseif(kernel == 2) % 2nd order
        if(solve_vel == 1)
            K{1} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{1}.*(1 - exp(-(fieldr_ext/epsilon).^2/2));
            K{1}(1,1) = 0;
            K{2} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{2}.*(1 - exp(-(fieldr_ext/epsilon).^2/2));
            K{2}(1,1) = 0;
        else
            G = -(log(fieldr_ext) + 1/2*expint((fieldr_ext/epsilon).^2/2))/(2*pi);
            G(1,1) = (gamma/2-log(sqrt(2)*epsilon))/(2*pi);
        end
    elseif(kernel == 4) % 4th order
        if(solve_vel == 1)
            K{1} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{1}.*(1 - (1 - (fieldr_ext/epsilon).^2/2).*exp(-(fieldr_ext/epsilon).^2/2));
            K{1}(1,1) = 0;
            K{2} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{2}.*(1 - (1 - (fieldr_ext/epsilon).^2/2).*exp(-(fieldr_ext/epsilon).^2/2));
            K{2}(1,1) = 0;
        else
            G = -(log(fieldr_ext) + 1/2*expint((fieldr_ext/epsilon).^2/2) - 1/2.*exp(-(fieldr_ext/epsilon).^2/2))/(2*pi);
            G(1,1) = (gamma/2 - log(sqrt(2)*epsilon) + 1/2)/(2*pi);
        end
    elseif(kernel == 6)% 6th order
        if(solve_vel == 1)
            K{1} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{1}.*(1 - (1 - (fieldr_ext/epsilon).^2 + 1/8*(fieldr_ext/epsilon).^4).*exp(-(fieldr_ext/epsilon).^2/2));
            K{1}(1,1) = 0;
            K{2} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{2}.*(1 - (1 - (fieldr_ext/epsilon).^2 + 1/8*(fieldr_ext/epsilon).^4).*exp(-(fieldr_ext/epsilon).^2/2));
            K{2}(1,1) = 0;
        else
            G = -(log(fieldr_ext) + 0.5*expint((fieldr_ext/epsilon).^2/2) - (3/4 - 1/8*(fieldr_ext/epsilon).^2).*exp(-(fieldr_ext/epsilon).^2/2))/(2*pi);
            G(1,1) = (gamma/2-log(sqrt(2)*epsilon) + 3/4)/(2*pi);
        end
    elseif(kernel == 8) % 8th order
        if(solve_vel == 1)
            K{1} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{1}.*(1 - (1 - 3/2*(fieldr_ext/epsilon).^2 + 3/8*(fieldr_ext/epsilon).^4 - 1/48*(fieldr_ext/epsilon).^6).*exp(-(fieldr_ext/epsilon).^2/2));
            K{1}(1,1) = 0;
            K{2} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{2}.*(1 - (1 - 3/2*(fieldr_ext/epsilon).^2 + 3/8*(fieldr_ext/epsilon).^4 - 1/48*(fieldr_ext/epsilon).^6).*exp(-(fieldr_ext/epsilon).^2/2));
            K{2}(1,1) = 0;
        else
            G = -(log(fieldr_ext) + 0.5*expint((fieldr_ext/epsilon).^2/2) - (11/12 - 7/24*(fieldr_ext/epsilon).^2 + 1/48*(fieldr_ext/epsilon).^4).*exp(-(fieldr_ext/epsilon).^2/2))/(2*pi);
            G(1,1) = (gamma/2-log(sqrt(2)*epsilon) + 11/12)/(2*pi);
        end
    elseif(kernel == 10) % 10th order
        if(solve_vel == 1)
            K{1} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{1}.*(1 - (1 - 2*(fieldr_ext/epsilon).^2 + 3/4*(fieldr_ext/epsilon).^4 - 1/12*(fieldr_ext/epsilon).^6 + 1/384*(fieldr_ext/epsilon).^8).*exp(-(fieldr_ext/epsilon).^2/2));
            K{1}(1,1) = 0;
            K{2} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{2}.*(1 - (1 - 2*(fieldr_ext/epsilon).^2 + 3/4*(fieldr_ext/epsilon).^4 - 1/12*(fieldr_ext/epsilon).^6 + 1/384*(fieldr_ext/epsilon).^8).*exp(-(fieldr_ext/epsilon).^2/2));
            K{2}(1,1) = 0;
        else
            G = -(log(fieldr_ext) + 0.5*expint((fieldr_ext/epsilon).^2/2) - (25/24 - 23/48*(fieldr_ext/epsilon).^2 + 13/192*(fieldr_ext/epsilon).^4 - 1/384*(fieldr_ext/epsilon).^6).*exp(-(fieldr_ext/epsilon).^2/2))/(2*pi);
            G(1,1) = (gamma/2-log(sqrt(2)*epsilon) + 25/24)/(2*pi);
        end
    elseif(kernel == 12) % 12th order
        if(solve_vel == 1)
            K{1} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{1}.*(1 - (1 - 5/2*(fieldr_ext/epsilon).^2 + 5/4*(fieldr_ext/epsilon).^4 - 5/24*(fieldr_ext/epsilon).^6 + 5/384*(fieldr_ext/epsilon).^8 - 1/3840*(fieldr_ext/epsilon).^10).*exp(-(fieldr_ext/epsilon).^2/2));
            K{1}(1,1) = 0;
            K{2} = -1./(2*pi*fieldr_ext.^2).*fieldx_ext{2}.*(1 - (1 - 5/2*(fieldr_ext/epsilon).^2 + 5/4*(fieldr_ext/epsilon).^4 - 5/24*(fieldr_ext/epsilon).^6 + 5/384*(fieldr_ext/epsilon).^8 - 1/3840*(fieldr_ext/epsilon).^10).*exp(-(fieldr_ext/epsilon).^2/2));
            K{2}(1,1) = 0;
        else
            G = -(log(fieldr_ext) + 0.5*expint((fieldr_ext/epsilon).^2/2) - (137/120 - 163/240*(fieldr_ext/epsilon).^2 + 137/960*(fieldr_ext/epsilon).^4 - 7/640*(fieldr_ext/epsilon).^6 + 1/3840*(fieldr_ext/epsilon).^8).*exp(-(fieldr_ext/epsilon).^2/2))/(2*pi);
            G(1,1) = (gamma/2-log(sqrt(2)*epsilon) + 137/120)/(2*pi);
        end
    elseif(kernel == 100) % Bessel
        if(solve_vel == 1)
            rho = pi*fieldr_ext/max(dx);
            K{1} = -(1 - besselj(0,rho))./(2*pi*fieldr_ext.^2).*fieldx_ext{1};
            K{1}(1,1) = 0;
            K{2} = -(1 - besselj(0,rho))./(2*pi*fieldr_ext.^2).*fieldx_ext{2};
            K{2}(1,1) = 0;                     
        else
            rho = pi*fieldr_ext/max(dx);
            G = 1/(2*pi)*(gamma - 1/8*(fieldr_ext/epsilon).^2*hypergeom([1,1],[2,2,2],- rho.^2/4) - log(2*rho) + log(rho.^2)/2);
            G(1,1) = -0.01845107351;
        end      
    else
        error('Specified kernel is not implemented')
    end   
    
    % Fourier transform integration kernel (normalised)
    if(solve_vel == 1)
        K_fft{1} = fft2(K{1}*dx(1)*dx(2));
        K_fft{2} = fft2(K{2}*dx(1)*dx(2));
    else
        G_fft = fft2(G*dx(1)*dx(2));
    end

%     clear G K (fieldr_ext/epsilon)
    %----------------------------------------------------------------------
    % Convolve and calculate the curl in Fourier space
    %----------------------------------------------------------------------
    if(solve_vel == 1)
        vel_fft{1} = vort_fft.*K_fft{2};
        vel_fft{2} = - vort_fft.*K_fft{1};
    else
        stream_fft = vort_fft.*G_fft;
        if(solve_vel == 2)
            vel_fft{1} = stream_fft.*(sqrt(-1)*fieldk{2});
            vel_fft{2} = -stream_fft.*(sqrt(-1)*fieldk{1});
        end
    end
    
    clear vort_fft K_fft G_fft fieldk
    %----------------------------------------------------------------------
    % Transform to physical space
    %----------------------------------------------------------------------
    if(solve_vel == 0)
        stream_ext = real(ifft2(stream_fft));
    else
        vel_ext{1} = real(ifft2(vel_fft{1}));
        vel_ext{2} = real(ifft2(vel_fft{2}));
    end
    
    %----------------------------------------------------------------------
    % Extract numerical solution
    %----------------------------------------------------------------------
    if(solve_vel == 0)
        stream = stream_ext(1:NX(1),1:NX(2));
        clear stream_fft stream_ext

        % Calculate diffference from analytical solution
        diff = sqrt((stream - stream_ana).^2) ...
            /max(max(abs(stream_ana)));
            
        roundoff = max(max(abs(stream_ana)))*eps;
    else
        vel{1} = vel_ext{1}(1:NX(1),1:NX(2));
        vel{2} = vel_ext{2}(1:NX(1),1:NX(2));
        
        clear vel_fft U2 V2
        
        % Calculate diffference from analytical solution
        diff = sqrt((vel{1} - vel_ana{1}).^2 + (vel{2} ...
            - vel_ana{2}).^2) ...
            /max(max(sqrt(vel_ana{1}.^2 + vel_ana{2}.^2)));
            
        roundoff = max(max(sqrt(vel_ana{1}.^2 + vel_ana{2}.^2)))*eps;
    end
    
    %----------------------------------------------------------------------
    % Error calculation (rms, max and L2-norm)
    %----------------------------------------------------------------------
    err_rms(N) = sqrt(sum(sum(diff.^2))/(NX(1)*NX(2)));
    err_max(N) = max(max(abs(diff)));
    err_L2(N)  = sum(sum(diff.^2))^(1/2)/(NX(1)*NX(2));
    
    dxs(N) = dx(1);
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
if(length(NXs) > 1)
    % Convergence test
    figure(10)
    set(gca,'FontSize',14)
    loglog(dxs,err_rms,colour,'LineWidth',2)
    hold on
     loglog(dxs,err_max,'g')
%     loglog(dxs,err_L2,'r')
    loglog(dxs,dxs.^2,'--k')
    loglog(dxs,dxs.^4,'--k')
    loglog(dxs,dxs.^6,'--k')
    loglog(dxs,dxs.^8,'--k')
    loglog(dxs,dxs.^10,'--k')
%     loglog(dxs,dxs.^22,'k')
    plot([min(dxs)/2 max(dxs)*2],[roundoff roundoff],'--r')

    xlim([min(dxs)/2 max(dxs)*2])
    title('Convergence test')
    xlabel('dx')
    ylabel('error')
%     legend('rms error','max error','L2 error','dx^2','dx^4','dx^6','dx^8','dx^{10}',1)
end
disp(['Round off: ' num2str(roundoff)]) 


% Initial and resulting field along (x,y) = (:,0)
figure(20)
set(gca,'FontSize',14)
plot(x{1},vort(:,ceil(NX(2)/2))/max(max(vort)),'b')
hold on
if(solve_vel == 0)
    plot(x{1},stream_ana(:,ceil(NX(2)/2))/max(max(stream_ana)),'r')
    plot(x{1},stream(:,ceil(NX(2)/2))/max(max(stream)),'--g')
    title('Initial vorticity and resulting stream field at (x,y) = (:,0)')
    xlabel('x')
    ylabel('normalised vort , stream')
%     legend('vorticity','analytic stream function','calculated stream function',1)
else
    plot(x{1},vel_ana{2}(:,ceil(NX(2)/2))/max(max(abs(vel_ana{2}))),'r')
    plot(x{1},vel{2}(:,ceil(NX(2)/2))/max(max(abs(vel_ana{2}))),'--g')
    title('Initial vorticity and resulting velocity field at (x,y) = (:,0)')
    xlabel('x')
    ylabel('normalised vort , u_2')
%     legend('vorticity','analytic velocity','calculated velocity',1)
end


% References
% [1] M. Perlman, On the accuracy of vortex methods, J. Comput. Phys. 59
%     (1985) 200-223
% [2] Beale, Majda (1985)
