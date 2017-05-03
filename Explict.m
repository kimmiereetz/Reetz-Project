clear all
clc

N=50;
dx=2*pi/(N+1);
dt=.99*dx^2/4;
t=0:dt:30;
k=dt/(dx^2);
c=1-4*k;
y=0:dx:2*pi;
x=y;
ny=length(y);
ts=length(t);
PHIab=cos(pi*y).*cosh(2*pi-y);
PSIab=y.^2.*sin(y/4);
U=[PHIab;zeros(N,N+2);PSIab];

[X,Y]=meshgrid(x,y);

%% 
% caluclate variables
Un=zeros(N+2,N+2);
for h=1:ts;
    for j=1:ny;
        for i=2:ny-1;
            if j==1
                Un(i,1)=k*U(i-1,1)+c*U(i,1)+k*U(i+1,1)+2*k*U(i,2);
            elseif j==ny;
                Un(i,j)=2*k*U(i,j-1)+k*U(i-1,j)+c*U(i,j)+k*U(i+1,j);
            else
                Un(i,j)=k*U(i,j-1)+k*U(i-1,j)+c*U(i,j)+k*U(i+1,j)+k*U(i,j+1);
            end
        end
        gi(h)=Un(N/2,N/2);
    end
    U=[U(1,:);Un(2:end-1,:);U(end,:)];
    
    if mod(h,5)==0;
        surf(X,Y,U)
        drawnow
        
    end
end
plot(t,gi)



