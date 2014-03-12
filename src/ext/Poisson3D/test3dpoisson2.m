clear
clc

%test function -laplace u=-2z(1-x+3y^2)/(1+x+y^2)^3 with exact solution
%u=z/(1+x+y^2)

err1=zeros(1,7);
err2=zeros(1,7);
h=zeros(1,7);
for m=1:7
nx=2^m;
ny=2^m;
nz=2^m;
hx=1/(nx+1);
hy=1/(ny+1);
hz=1/(nz+1);
x=linspace(0,1,nx+2);
y=linspace(0,1,ny+2);
z=linspace(0,1,nz+2);
f=zeros(nx+2,ny+2,nz+2);
uexact=zeros(nx+2,ny+2,nz+2);
for k=1:nz+2
    for j=1:ny+2
        for i=1:nx+2
            f(i,j,k)=-2*z(k)*(1/(1+x(i)+y(j)^2)^3-1/(1+x(i)+y(j)^2)^2+4*y(j)^2/(1+x(i)+y(j)^2)^3);
            uexact(i,j,k)=z(k)/(1+x(i)+y(j)^2);
        end
    end
end

u1=poisson3d(f,hx^2,hy^2,hz^2,uexact(1,:,:),uexact(nx+2,:,:),uexact(:,1,:),uexact(:,ny+2,:),uexact(:,:,1),uexact(:,:,nz+2),2);
u2=poisson3d(f,hx^2,hy^2,hz^2,uexact(1,:,:),uexact(nx+2,:,:),uexact(:,1,:),uexact(:,ny+2,:),uexact(:,:,1),uexact(:,:,nz+2),4);

err1(m)=max(max(max(abs(uexact(2:nx+1,2:ny+1,2:nz+1)-u1))));
err2(m)=max(max(max(abs(uexact(2:nx+1,2:ny+1,2:nz+1)-u2))));
h(m)=hx;
end


figure(1)
loglog(h,[err1;err2],'*-')
title('h vs err(\infty norm)')
xlabel('h')
ylabel('err')
legend('2nd-order scheme','4th-order scheme','Location','Best')

slope1=zeros(1,6);
slope2=zeros(1,6);
for m=2:7
    slope1(m-1)=log(err1(m)/err1(m-1))/log(h(m)/h(m-1));
    slope2(m-1)=log(err2(m)/err2(m-1))/log(h(m)/h(m-1));
end

figure(2)
plot(1:6,[slope1;slope2],'*-')
title('slope of loglog graph')
xlabel('step number')
ylabel('slope')
legend('2nd-order scheme','4th-order scheme','Location','Best')