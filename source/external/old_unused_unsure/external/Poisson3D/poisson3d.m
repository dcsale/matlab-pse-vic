function u=poisson3d(f,hx2,hy2,hz2,bx1,bx2,by1,by2,bz1,bz2,order)
% A 3D Fast Poisson Solver for the equation -Laplace(u)=f in a cuboid with the
% Dirichlet boundary condition.
%
% written by Weiguo Gao (wggao@fudan.edu.cn) and Yuan Cao
%
% poisson3d(f,hx2,hy2,hz2,bx1,bx2,by1,by2,bz1,bz2,order):
% - f is a 3D array of the values of the equation's right hand function on
%   mesh points.
% - hx2,hy2 and hz2 are the squares of the distance between two adjacent
%   points in x,y,z directions.
% - bx1,bx2,by1,by2,bz1,bz2 are the corresponding boundary conditions of the
%   equation.
% - order=2 or 4 is the parameter for users to choose a scheme. Set order=2
%   for a 2nd-order scheme, and set order=4 for a 4th order one.
%
% This file requires three additional files:
% my3dDST.m/my3diDST.m  -  calls discrete sin transform/ inverse
%                          discrete sin transform directly for 3D arrays.
% secdiff               -  forms the modified right-hand side of the problem
%                          when the 4th-order scheme is resorted to.

[nx ny nz]=size(f);
change=0;

if order==2
%deal with the boundary conditions.
f(2,:,:)=f(2,:,:)+bx1/hx2;
f(nx-1,:,:)=f(nx-1,:,:)+bx2/hx2;

f(:,2,:)=f(:,2,:)+by1/hy2;
f(:,ny-1,:)=f(:,ny-1,:)+by2/hy2;

f(:,:,2)=f(:,:,2)+bz1/hz2;
f(:,:,nz-1)=f(:,:,nz-1)+bz2/hz2;

u=f(2:nx-1,2:ny-1,2:nz-1);
clear bx1 bx2 by1 by2 bz1 bz2 f

nx=nx-2;
ny=ny-2;
nz=nz-2;

%Choose the longest direction to solve the tridiagonal linear equations to
%make the calculation faster.
if nx>ny && nx>nz
    change=1;
    changeindex=[3 2 1];
    f=permute(f,changeindex);
    [nx ny nz]=size(f);
    temp=hx2;
    hx2=hz2;
    hz2=temp;
elseif ny>nx && ny>nz
    change=1;
    changeindex=[1 3 2];
    f=permute(f,changeindex);
    [nx ny nz]=size(f);
    temp=hy2;
    hy2=hz2;
    hz2=temp;
end

%inverse discrete sine transform in the other two directions.
u=my3diDST(u);
u=permute(my3diDST(permute(u,[2 1 3])),[2 1 3]);

%solve the tridiagonal linear equations in the longest direction.
lambdax=4/hx2*(sin((1:nx)*pi/2/(nx+1))).^2;
lambday=4/hy2*(sin((1:ny)*pi/2/(ny+1))).^2;
trid=sparse([1:nz-1 2:nz 1:nz]',[2:nz 1:nz-1 1:nz]',-ones(3*nz-2,1)/hz2);
for i=1:nx
    for j=1:ny
        trid(1:nz+1:nz*nz)=lambdax(i)+lambday(j)+2/hz2;
        u(i,j,:)=trid\reshape(u(i,j,:),nz,1);
    end
end

%call discrete sine transform to get the solution.
u=my3dDST(u);
u=permute(my3dDST(permute(u,[2 1 3])),[2 1 3]);


elseif order==4
nx=nx-2;ny=ny-2;nz=nz-2;

%forms the modified right-hand side of the problem.
f=f(2:nx+1,2:ny+1,2:nz+1)-(secdiff(f)+permute(secdiff(permute(f,[2,1,3])),[2,1,3])+permute(secdiff(permute(f,[3,1,2])),[2,3,1]))/12;

%deal with the boundary conditions.
f(1,:,:)=f(1,:,:)-(-2/3/hx2+1/6/hy2+1/6/hz2)*bx1(1,2:ny+1,2:nz+1);
f(nx,:,:)=f(nx,:,:)-(-2/3/hx2+1/6/hy2+1/6/hz2)*bx2(1,2:ny+1,2:nz+1);
f(:,1,:)=f(:,1,:)-(-2/3/hy2+1/6/hx2+1/6/hz2)*by1(2:nx+1,1,2:nz+1);
f(:,ny,:)=f(:,ny,:)-(-2/3/hy2+1/6/hx2+1/6/hz2)*by2(2:nx+1,1,2:nz+1);
f(:,:,1)=f(:,:,1)-(-2/3/hz2+1/6/hx2+1/6/hy2)*bz1(2:nx+1,2:ny+1,1);
f(:,:,nz)=f(:,:,nz)-(-2/3/hz2+1/6/hx2+1/6/hy2)*bz2(2:nx+1,2:ny+1,1);

f(1,:,:)=f(1,:,:)+1/12*(1/hx2+1/hy2)*(bx1(1,1:ny,2:nz+1)+bx1(1,3:ny+2,2:nz+1))+1/12*(1/hx2+1/hz2)*(bx1(1,2:ny+1,1:nz)+bx1(1,2:ny+1,3:nz+2));
f(nx,:,:)=f(nx,:,:)+1/12*(1/hx2+1/hy2)*(bx2(1,1:ny,2:nz+1)+bx2(1,3:ny+2,2:nz+1))+1/12*(1/hx2+1/hz2)*(bx2(1,2:ny+1,1:nz)+bx2(1,2:ny+1,3:nz+2));

f(2:nx,1,:)=f(2:nx,1,:)+1/12*(1/hx2+1/hy2)*by1(2:nx,1,2:nz+1);
f(1:nx-1,1,:)=f(1:nx-1,1,:)+1/12*(1/hx2+1/hy2)*by1(3:nx+1,1,2:nz+1);
f(:,1,:)=f(:,1,:)+1/12*(1/hy2+1/hz2)*(by1(2:nx+1,1,1:nz)+by1(2:nx+1,1,3:nz+2));

f(2:nx,ny,:)=f(2:nx,ny,:)+1/12*(1/hx2+1/hy2)*by2(2:nx,1,2:nz+1);
f(1:nx-1,ny,:)=f(1:nx-1,ny,:)+1/12*(1/hx2+1/hy2)*by2(3:nx+1,1,2:nz+1);
f(:,ny,:)=f(:,ny,:)+1/12*(1/hy2+1/hz2)*(by2(2:nx+1,1,1:nz)+by2(2:nx+1,1,3:nz+2));

f(2:nx,:,1)=f(2:nx,:,1)+1/12*(1/hx2+1/hz2)*bz1(2:nx,2:ny+1,1);
f(1:nx-1,:,1)=f(1:nx-1,:,1)+1/12*(1/hx2+1/hz2)*bz1(3:nx+1,2:ny+1,1);
f(:,2:ny,1)=f(:,2:ny,1)+1/12*(1/hy2+1/hz2)*bz1(2:nx+1,2:ny,1);
f(:,1:ny-1,1)=f(:,1:ny-1,1)+1/12*(1/hy2+1/hz2)*bz1(2:nx+1,3:ny+1,1);

f(2:nx,:,nz)=f(2:nx,:,nz)+1/12*(1/hx2+1/hz2)*bz2(2:nx,2:ny+1,1);
f(1:nx-1,:,nz)=f(1:nx-1,:,nz)+1/12*(1/hx2+1/hz2)*bz2(3:nx+1,2:ny+1,1);
f(:,2:ny,nz)=f(:,2:ny,nz)+1/12*(1/hy2+1/hz2)*bz2(2:nx+1,2:ny,1);
f(:,1:ny-1,nz)=f(:,1:ny-1,nz)+1/12*(1/hy2+1/hz2)*bz2(2:nx+1,3:ny+1,1);
clear bx1 bx2 by1 by2 bz1 bz2

%Choose the longest direction to solve the tridiagonal linear equations to
%make the calculation faster.
if nx>ny && nx>nz
    change=1;
    changeindex=[3 2 1];
    f=permute(f,changeindex);
    [nx ny nz]=size(f);
    temp=hx2;
    hx2=hz2;
    hz2=temp;
elseif ny>nx && ny>nz
    change=1;
    changeindex=[1 3 2];
    f=permute(f,changeindex);
    [nx ny nz]=size(f);
    temp=hy2;
    hy2=hz2;
    hz2=temp;
end

%inverse discrete sine transform in the other two directions.
u=my3diDST(f);
clear f
u=permute(my3diDST(permute(u,[2 1 3])),[2 1 3]);

%solve the tridiagonal linear equations in the longest direction.
lambdax=4/hx2*(sin((1:nx)*pi/2/(nx+1))).^2;
lambday=4/hy2*(sin((1:ny)*pi/2/(ny+1))).^2;
trid=sparse([1:nz-1 2:nz 1:nz]',[2:nz 1:nz-1 1:nz]',-ones(3*nz-2,1)/hz2);
for i=1:nx
    for j=1:ny
        p=1-(hx2+hz2)/12*lambdax(i)-(hy2+hz2)/12*lambday(j);
        trid(1:nz+1:nz*nz)=2/hz2+(lambdax(i)+lambday(j)-(hx2+hy2)/12*lambdax(i)*lambday(j))/p;
        u(i,j,:)=trid\reshape(u(i,j,:),nz,1)/p;
    end
end

%call discrete sine transform to get the solution.
u=my3dDST(u);
u=permute(my3dDST(permute(u,[2 1 3])),[2 1 3]);
end

if change==1
    u=permute(u,changeindex);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=my3diDST(x)
%Discrete sin transform for 3D arrays. (c.f. dst)
y=2/(size(x,1)+1)*my3dDST(x);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=my3dDST(x)
%Discrete sin transform for 3D arrays. (c.f. dst)
n=size(x,1);
sizex=size(x);
m=length(sizex);
y=zeros(sizex);
switch m
    case 3
        xx=zeros([2*(n+1),sizex(2:3)]);
        xx(2:n+1,:,:)=x;
        xx(n+3:2*(n+1),:,:)=-x(n:-1:1,:,:);
        xx=imag(fft(xx));
        y=xx(2:n+1,:,:)/(-2);
    case 2
        xx=zeros(2*(n+1),sizex(2));
        xx(2:n+1,:)=x;
        xx(n+3:2*(n+1),:)=-x(n:-1:1,:);
        xx=imag(fft(xx));
        y=xx(2:n+1,:)/(-2);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ff=secdiff(f)
%forms the modified right-hand side of the problem when the 4th-order
%scheme is resorted to.
[nx ny nz]=size(f);
nx=nx-2;ny=ny-2;nz=nz-2;
T=sparse([1:nx 1:nx 1:nx]',[1:nx 2:nx+1 3:nx+2]',[-ones(nx,1);2*ones(nx,1);-ones(nx,1)]);
ff=zeros(nx,ny,nz);
for k=1:nz
    ff(:,:,k)=T*f(:,2:ny+1,k+1);
end

return