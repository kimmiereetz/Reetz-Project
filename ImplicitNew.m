clear all
clc

N=40;
dx=2*pi/(N-1);
dt=.9*dx^2/4;
x=0:dx:2*pi;
y=x;
t=0:dt:30;
nx=length(x);
ts=length(t);
%%

a=dt/(2*dx^2);
b=(1+2*a);
c=(1-2*a);
[X,Y]=meshgrid(x,y);
%%

PHIab=cos(pi*y).*cosh(2*pi-y);
PSIab=y.^2.*sin(y/4);
U=[PHIab;zeros(nx-2,nx);PSIab];

%%
% Tridiag: Tridiagonal equation solver banded system
%   x = Tridiag(e,f,g,r): Tridiagonal system solver. 
% input: 
%   e = subdiagonal vector 
%   f = diagonal vector 
%   g = superdiagonal vector 
%   r = right hand side vector 
% output: 
%   x = solution vector
e=[-a*ones(1,N-2) 0];
f=[ 1 b*ones(1, N-2) 1];
g=[ 0 -a*ones(1, N-2)];
r=zeros(1,N);

e2=[-a*ones(1,N-2) -2*a];
f2=b*ones(N);
g2=[-2*a -a*ones(1,N-2)];
r2=zeros(1,N);

for h=1:ts;
    for j=1:nx;
        for i=2:nx-1;
            if j==1
                r(i)=c*U(i,1)+2*a*U(i,2);
            elseif j==nx
                r(i)=2*a*U(i,nx-1)+c*U(i,nx);
            else
                r(i)=a*U(i,j-1)+c*U(i,j)+a*U(i,j+1);
            end
        end
        r=[U(1,j) r(2:nx-1) U(nx,j)];
        x=tridiag(e,f,g,r);
        U(:,j)=x;
    end
    


%%


    for i=2:nx-1;
        for j=1:nx;
             r2(j)=a*U(i-1,j)+c*U(i,j)+a*U(i+1,j);
        end
        x2=tridiag(e2,f2,g2,r2);
        U(i,:)=x2;
    end
    center(h)=U(N/2,N/2);
        figure (1)
    if mod(h,5)==0;
        surf(X,Y,U)
        drawnow
    end
        
end
figure (2)
plot(t,center)
gi=norm(U);
norm(U);