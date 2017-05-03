clear all
clc

N=20;
dx=2*pi/(N-1);
dt=.9*dx^2/4;
x=0:dx:2*pi;
y=x;
t=0:dt:10;
nx=length(x);
%%

a=dt/(2*dx^2);
b=(1+2*a);
c=(1-2*a);

%%

PHIab=cos(pi*y).*cosh(2*pi-y);
PSIab=y.^2.*sin(y/4);
U=[PHIab;zeros(nx-2,nx);PSIab];

%%
r=zeros(1,N);
for j=1:nx;
    for i=1:nx;
        if j==1
            r(i)=c*U(i,1)+2*a*U(i,2);
        elseif j==nx
            r(i)=2*a*U(i,nx-1)+c*U(i,nx);
        else
            r(i)=a*U(i,j-1)+c*U(i,j)+a*U(i,j+1);
        end
    end
end



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


x=tridiag(e,f,g,r);