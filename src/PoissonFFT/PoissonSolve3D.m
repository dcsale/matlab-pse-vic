% function [field_out] = PoissonSolve3D(vort, NX, Mesh.k, fieldxf_ext, kernel, solve_vel, runTests)
function field_out = PoissonSolve3D(vort, Mesh, Sim)
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
%* Changelog:    2/22/2014 - Danny Clay Sale (dsale@uw.edu)
%                * split original scripts into callable functions
%*               * Added debugging flags for plots and tests                          
%*               * GPU acceleration of ndgrid() and fftn()
%*               * renamed some variables to clarify field and coordinates
%*               
%**************************************************************************

do_timing = true;
if do_timing
    tic
end

%----------------------------------------------------------------------
% Fourier transform vorticity field
%----------------------------------------------------------------------
vort_ext = cell(1, 3);
vort_fft = cell(1, 3);
% for i = 1:3
parfor i = 1:3    
    % the extended vorticity field
    vort_ext{i}                                         = zeros(2*Mesh.NX(1), 2*Mesh.NX(2), 2*Mesh.NX(3));
    vort_ext{i}(1:Mesh.NX(1),1:Mesh.NX(2),1:Mesh.NX(3)) = vort{i};
    vort_fft{i}                                         = fftn(vort_ext{i});
end
clear vort_ext

%----------------------------------------------------------------------
% Integration kernel
%----------------------------------------------------------------------
epsilon = Sim.alpha*max(Mesh.dx); % smoothing radius
rho     = Mesh.rf_ext/epsilon; % normalised radial distance
if(Sim.kernel == 0) % 0th order (No regularisation)
    if(Sim.solve_vel == 1)
        K{1} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{1};
        K{1}(1,1,1) = 0;
        K{2} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{2};
        K{2}(1,1,1) = 0;
        K{3} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{3};
        K{3}(1,1,1) = 0;
    else
        G = 1./(4*pi*Mesh.rf_ext);
        G(1,1,1) = 1;
    end
elseif(Sim.kernel == 2) % 2nd order
    if(Sim.solve_vel == 1)
        K{1} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{1}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{1}(1,1,1) = 0;
        K{2} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{2};
        K{2}(1,1,1) = 0;
        K{3} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{3};
        K{3}(1,1,1) = 0;            
    else
        G = 1./(4*pi*Mesh.rf_ext).*erf(rho/sqrt(2));
        G(1,1,1) = sqrt(2)/(4*pi^(3/2)*epsilon);
    end
elseif(Sim.kernel == 4) % 4th order
    if(Sim.solve_vel == 1)
        K{1} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{1}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{1}(1,1,1) = 0;
        K{2} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{2}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{2}(1,1,1) = 0;
        K{3} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{3}.*((-2*rho + rho.^3).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{3}(1,1,1) = 0;
    else
        G = 1./(4*pi*Mesh.rf_ext).*((rho).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
        G(1,1,1) = 3*sqrt(2)/(8*pi^(3/2)*epsilon);
    end
elseif(Sim.kernel == 6) % 6th order
    if(Sim.solve_vel == 1)
        K{1} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{1}.*((-2*rho + 9/4*rho.^3 - 1/4*rho.^5).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{1}(1,1,1) = 0;
        K{2} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{2}.*((-2*rho + 9/4*rho.^3 - 1/4*rho.^5).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{2}(1,1,1) = 0;
        K{3} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{3}.*((-2*rho + 9/4*rho.^3 - 1/4*rho.^5).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{3}(1,1,1) = 0;            
    else
        G = 1./(4*pi*Mesh.rf_ext).*((7/4*rho - 1/4*rho.^3).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
        G(1,1,1) = 15*sqrt(2)/(32*pi^(3/2)*epsilon);
    end
elseif(Sim.kernel == 8) % 8th order
    if(Sim.solve_vel == 1)
        K{1} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{1}.*((-2*rho + 89/24*rho.^3 - 20/24*rho.^5 + 1/24*rho.^7).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{1}(1,1,1) = 0;
        K{2} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{2}.*((-2*rho + 89/24*rho.^3 - 20/24*rho.^5 + 1/24*rho.^7).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{2}(1,1,1) = 0;
        K{3} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{3}.*((-2*rho + 89/24*rho.^3 - 20/24*rho.^5 + 1/24*rho.^7).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{3}(1,1,1) = 0;            
    else
        G = 1./(4*pi*Mesh.rf_ext).*((19/8*rho - 2/3*rho.^3 + 1/24*rho.^5).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
        G(1,1,1) = 35*sqrt(2)/(64*pi^(3/2)*epsilon);
    end
elseif(Sim.kernel == 10) % 10th order
    if(Sim.solve_vel == 1)
        K{1} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{1}.*((-2*rho + 1027/192*rho.^3 - 349/192*rho.^5 + 35/192*rho.^7 - 1/192*rho.^9).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{1}(1,1,1) = 0;
        K{2} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{2}.*((-2*rho + 1027/192*rho.^3 - 349/192*rho.^5 + 35/192*rho.^7 - 1/192*rho.^9).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{2}(1,1,1) = 0;
        K{3} = - 1./(4*pi*Mesh.rf_ext.^3).*Mesh.xf_ext{3}.*((-2*rho + 1027/192*rho.^3 - 349/192*rho.^5 + 35/192*rho.^7 - 1/192*rho.^9).*exp(-1/2*rho.^2)/sqrt(2*pi) + erf(rho/sqrt(2)));
        K{3}(1,1,1) = 0;            
    else
        G = 1./(4*pi*Mesh.rf_ext).*((187/64*rho - 233/192*rho.^3 + 29/192*rho.^5 - 1/192*rho.^7).*exp(-1/2*rho.^2)/(sqrt(2*pi)) + erf(rho/sqrt(2)));
        G(1,1,1) = 315*sqrt(2)/(512*pi^(3/2)*epsilon);
    end      
else
    error('Specified Sim.kernel is not implemented')
end

% Fourier transform integration kernel (normalised)
if(Sim.solve_vel == 1)
    K_fft{1} = fftn(K{1}*Mesh.dx(1)*Mesh.dx(2)*Mesh.dx(3));
    K_fft{2} = fftn(K{2}*Mesh.dx(1)*Mesh.dx(2)*Mesh.dx(3));
    K_fft{3} = fftn(K{3}*Mesh.dx(1)*Mesh.dx(2)*Mesh.dx(3));
else
    G_fft = fftn(G*Mesh.dx(1)*Mesh.dx(2)*Mesh.dx(3));
end

%----------------------------------------------------------------------
% Convolve and calculate the curl in Fourier space
%----------------------------------------------------------------------
if(Sim.solve_vel == 1)
    vel_fft{1} = vort_fft{3}.*K_fft{2} - vort_fft{2}.*K_fft{3};
    vel_fft{2} = vort_fft{1}.*K_fft{3} - vort_fft{3}.*K_fft{1};
    vel_fft{3} = vort_fft{2}.*K_fft{1} - vort_fft{1}.*K_fft{2};
else
    stream_fft{1} = vort_fft{1}.*G_fft;
    stream_fft{2} = vort_fft{2}.*G_fft;
    stream_fft{3} = vort_fft{3}.*G_fft;
    if(Sim.solve_vel == 2)
        vel_fft{1} = stream_fft{3}.*(sqrt(-1)*Mesh.kf{2}) - stream_fft{2}.*(sqrt(-1)*Mesh.kf{3});
        vel_fft{2} = stream_fft{1}.*(sqrt(-1)*Mesh.kf{3}) - stream_fft{3}.*(sqrt(-1)*Mesh.kf{1});
        vel_fft{3} = stream_fft{2}.*(sqrt(-1)*Mesh.kf{1}) - stream_fft{1}.*(sqrt(-1)*Mesh.kf{2});
    end
end

%----------------------------------------------------------------------
% Transform back to physical space
%----------------------------------------------------------------------
if(Sim.solve_vel == 0)
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
if Sim.solve_vel == 0 
    stream{1} = stream_ext{1}(1:Mesh.NX(1), 1:Mesh.NX(2), 1:Mesh.NX(3));
    stream{2} = stream_ext{2}(1:Mesh.NX(1), 1:Mesh.NX(2), 1:Mesh.NX(3));
    stream{3} = stream_ext{3}(1:Mesh.NX(1), 1:Mesh.NX(2), 1:Mesh.NX(3));  
    field_out = stream;
else
    vel{1}    = vel_ext{1}(1:Mesh.NX(1), 1:Mesh.NX(2), 1:Mesh.NX(3));
    vel{2}    = vel_ext{2}(1:Mesh.NX(1), 1:Mesh.NX(2), 1:Mesh.NX(3));
    vel{3}    = vel_ext{3}(1:Mesh.NX(1), 1:Mesh.NX(2), 1:Mesh.NX(3));
    field_out = vel;
end

if do_timing
%     user = memory;
%     mem  = user.MemUsedMATLAB;
    Mesh.NX 
    toc 
%     mem/1e9
%     memory
end

end % function PoissonSolve3D

% function Results = testSolver()
%     
% %% NOTE: put Mads original way of calling the solver in here (for convergence studies)
% 
% ----------------------------------------------------------------------
% UNIT TESTS !!! 
% run the original test examples, analytical solutions & plots
% ----------------------------------------------------------------------
% if runTests(Sim.solve_vel)
%     
%     % Calculate difference from analytical solution
%     if Sim.solve_vel == 0
%         diff = sqrt((stream{1} - stream_ana{1}).^2 + (stream{2} ...
%             - stream_ana{2}).^2 + (stream{3} - stream_ana{3}).^2) ...
%             /max(max(max(sqrt(stream_ana{1}.^2 + stream_ana{2}.^2 + stream_ana{3}.^2))));
%     else
%         diff = sqrt((vel{1} - vel_ana{1}).^2 + (vel{2} ...
%             - vel_ana{2}).^2 + (vel{3} - vel_ana{3}).^2) ...
%             /max(max(max(sqrt(vel_ana{1}.^2 + vel_ana{2}.^2 + vel_ana{3}.^2))));
%     end
% 
%     %----------------------------------------------------------------------
%     % Error calculation (rms, max and L2-norm)
%     %----------------------------------------------------------------------
%     err_rms(N) = sqrt(sum(sum(sum(diff.^2)))/(Mesh.NX(1)*Mesh.NX(2)*Mesh.NX(3)));
%     err_max(N) = max(max(max(abs(diff))));    
%     err_L2(N)  = sum(sum(sum(diff.^2)))^(1/2)/(Mesh.NX(1)*Mesh.NX(2)*Mesh.NX(3));
% 
%     dxs(N)     = dx(1);
%     %--------------------------------------------------------------------------
%     % Plot results
%     %--------------------------------------------------------------------------
%     if(length(Mesh.NXs) > 1)
%         % Convergence test
%         figure(10)
%         set(gca,'FontSize',14)
%         loglog(dxs,err_L2,'b')
%         hold on
%         loglog(dxs,err_max,'g')
%         loglog(dxs,err_rms,'r')
%         loglog(dxs,dxs.^2,'k')
%         loglog(dxs,dxs.^4,'k')
%         loglog(dxs,dxs.^6,'k')
%         loglog(dxs,dxs.^8,'k')
%         loglog(dxs,dxs.^10,'k')
%         xlim([min(dxs)/2 max(dxs)*2])
%         title('Convergence test')
%         xlabel('dx')
%         ylabel('error')
%         legend('rms error','max error','L2 error','dx^2','dx^4','dx^6','dx^8','dx^{10}',1)
%     end
% 
% 
%     %----------------------------------------------------------------------
%     % Analytical solution
%     %----------------------------------------------------------------------
%     if(Sim.solve_vel == 0)
%         disp(['Solving stream function using kernel G' num2str(Sim.kernel) ' at Mesh.NX = ' num2str(Mesh.NX)])
%         stream_ana{1} = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
%         stream_ana{2} = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
%         stream_ana{3} = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
% 
%         if(testcase == 1)
%             % Bump function: SPHERICAL SCALAR FIELD
%             stream_ana{2}(fieldr_cen < R) = exp(-c./(1 - fieldr_cen(fieldr_cen < R).^2/R^2));
% 
%         elseif(testcase == 2)
%             % Bump function: VORTEX RING (in the xy-plane)
%             stream_ana_mag = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
%             stream_ana_mag(fieldphi_cen < R) = exp(-c./(1 - (fieldr_cen(fieldphi_cen < R).^2 ...
%                 + R^2 + fieldx_cen{3}(fieldphi_cen < R).^2 - 2*R*fieldr_cen(fieldphi_cen < R))/R^2));
%             stream_ana{1} = -sin(fieldtheta).*stream_ana_mag;
%             stream_ana{2} = cos(fieldtheta).*stream_ana_mag;
%             stream_ana{3} = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
%         end
%     else
%         if(Sim.solve_vel == 1)
%             disp(['Solving velocity using kernel K' num2str(Sim.kernel) ' at Mesh.NX = ' num2str(Mesh.NX)])
%         else
%             disp(['Solving velocity using kernel G' num2str(Sim.kernel) ' at Mesh.NX = ' num2str(Mesh.NX)])
%         end
%         vel_ana{1} = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
%         vel_ana{2} = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
%         vel_ana{3} = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
% 
%         if(testcase == 1)
%             % Bump function: SPHERICAL SCALAR FIELD
%             error('Velocity field is not available for this test function - set Sim.solve_vel = 0')
%         elseif(testcase == 2)
%             % Bump function VORTEX RING (in the xy-plane)
%             vel_ana_r = zeros(Mesh.NX(1),Mesh.NX(2),Mesh.NX(3));
%             vel_ana_r(fieldphi_cen < R) = 2*c*R^2*fieldx_cen{3}(fieldphi_cen < R) ...
%                 .* exp(-c*R^2./(2*R*fieldr_cen(fieldphi_cen < R) - ...
%                 fieldr_cen(fieldphi_cen < R).^2 - fieldx_cen{3}(fieldphi_cen < R).^2)) ...
%                 .*(2*R*fieldr_cen(fieldphi_cen < R) - fieldr_cen(fieldphi_cen < R).^2 ...
%                 - fieldx_cen{3}(fieldphi_cen < R).^2).^(-2);
% 
%             vel_ana{1} = cos(fieldtheta).*vel_ana_r;
%             vel_ana{2} = sin(fieldtheta).*vel_ana_r;
% 
%             vel_ana{3}(fieldphi_cen < R) = exp(-c*R^2./(2*R*fieldr_cen(fieldphi_cen < R) - ...
%                 fieldr_cen(fieldphi_cen < R).^2 - fieldx_cen{3}(fieldphi_cen < R).^2)) ...
%                 .*(4*R^2*fieldr_cen(fieldphi_cen < R).^2 ...
%                 - 4*R*fieldr_cen(fieldphi_cen < R).^3 ...
%                 - 4*R*fieldr_cen(fieldphi_cen < R).*fieldx_cen{3}(fieldphi_cen < R).^2 ...
%                 + fieldr_cen(fieldphi_cen < R).^4 ...
%                 + 2*fieldr_cen(fieldphi_cen < R).^2.*fieldx_cen{3}(fieldphi_cen < R).^2 ...
%                 + fieldx_cen{3}(fieldphi_cen < R).^4 ...
%                 + 2*c*R^3*fieldr_cen(fieldphi_cen < R) ...
%                 - 2*c*R^2*fieldr_cen(fieldphi_cen < R).^2) ...
%                 .*(2*R*fieldr_cen(fieldphi_cen < R) - fieldr_cen(fieldphi_cen < R).^2 ...
%                 - fieldx_cen{3}(fieldphi_cen < R).^2).^(-2).*fieldr_cen(fieldphi_cen < R).^(-1);
%         end        
%     end
% 
%     clear fieldx_cen{1} fieldr_cen
%     
%     if plot_tests
%         
%         % Initial and resulting field along (x,y,z) = (:,0,0)
%         figure(20)
%         set(gca,'FontSize',14)
%         plot(x{1},vort{2}(:,Mesh.NX(2)/2,Mesh.NX(3)/2)/max(max(max(vort{2}))),'b')
%         hold on
%         if(Sim.solve_vel == 0)
%             plot(x{1},stream_ana{2}(:,Mesh.NX(2)/2,Mesh.NX(3)/2)/max(max(max(stream_ana{2}))),'r')
%             plot(x{1},stream{2}(:,Mesh.NX(2)/2,Mesh.NX(3)/2)/max(max(max(stream{2}))),'--g')
%             title('Initial vorticity and resulting stream field at (x,y,z) = (:,0,0)')
%             xlabel('x')
%             ylabel('normalised \omega_2 , \psi_3')
%             legend('vorticity','analytic stream function','calculated stream function',4)
%         else
%             plot(x{1},vel_ana{3}(:,Mesh.NX(2)/2,Mesh.NX(3)/2)/max(max(max(vel_ana{3}))),'r')
%             plot(x{1},vel{3}(:,Mesh.NX(2)/2,Mesh.NX(3)/2)/max(max(max(vel{3}))),'--g')  
%             title('Initial vorticity and resulting velocity field at (x,y,z) = (:,0,0)')
%             xlabel('x')
%             ylabel('normalised \omega_2 , u_3')
%             legend('vorticity','analytic velocity','calculated velocity',4)
%         end
%         
%     end
%     
% end
% 
% end % function testSolver

