clear all
clc

N=40;
    
dx=2*pi/(N-1);

dt=.9*dx^2/4;
t=0:dt:30;
y=0:dx:2*pi;
x=y;

k=dt/(dx^2);
c=1-4*k;

ny=length(y);
ts=length(t);
PHIab=cos(pi*y).*cosh(2*pi-y);
PSIab=y.^2.*sin(y/4);
% Creating the intial U grid
% Dirichlet boundary condition u(x=a_x,y)=?_ab(y)=cos[?*y].*cosh(2?-y)
% is imposed at U(1,:) 
% Dirichlet boundary condition u(x=b_x,y)=?_ab(y)=?y?.^2 sin(y/4) is
% imposed at U(ny,:)
% The intial condition u(x,y,t=0)= u_0 (x,y)=0 is imposed at all other
% points in the intial grid
U=[PHIab;zeros(N-2,N);PSIab];
% Allocating memory for variables
Un=zeros(ny,ny);
center=zeros(1,ts);
[X,Y]=meshgrid(x,y);


%% 
% caluclate variables


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
            center(h)=Un(N/2,N/2);
        end
        U=[U(1,:); Un(2:end-1,:) ; U(end,:)];
    
        
         if mod(h,10)==0;
             figure (1)
             surf(X,Y,U)
             drawnow
             
         end
    end

figure (2)
plot(t,center)
%gi=norm(U,2);


display('Done')