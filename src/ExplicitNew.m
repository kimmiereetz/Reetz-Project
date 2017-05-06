clear all
clc
% Setting up grid number of nodes N determines the fine-ness of the mesh
% the larger N the finer the mesh

N=40;
% Based on the number of nodes the increments for dx dy and dt are
% determined    
dx=2*pi/(N-1);
% Stability condition enforced
dt=.9*dx^2/4;
% Based on the increment size for the time and space dimensions and their
% intervals the vectors for x, y, and time are created
t=0:dt:30;
y=0:dx:2*pi;
x=y;
% The constants determined from discritization are calculated befoe the
% loop to optimize code

k=dt/(dx^2);
c=1-4*k;
% For ease of operating loops the number of elements in each x y and time
% vector are calculated
nx=length(x);
ts=length(t);

ny=length(y);
% The Dirichlet boundary conditions were reduced and are evaluated to 
% produce vectors.

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
% This mesh is required for plotting the surface solution 
[X,Y]=meshgrid(x,y);


%% 

% The loops for solving the problem are ready to run. Time is the outer
% most loop folloeing by the y direction and x direction

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
% The centrally located grid node is tracked for use of steady state
% determination and diffusion confirmation
           
            center(h)=Un(N/2,N/2);
        end
        U=[U(1,:); Un(2:end-1,:) ; U(end,:)];
    
% Throughout the time steps plots are generated to confirm diffusion
       
         if mod(h,10)==0;
             figure (1)
             surf(X,Y,U)
             drawnow
             
         end
    end

figure (2)
plot(t,center)
%gi=norm(U,2);% The L2 Norm of the final U matrix is used for grid convergence. The L2
% Norms at Nodes of N 2N 4N 16N are calculted until the solution truncates 
% The minimum difference in L2 Norms is calculated to determine the sweet
% spot for the numerical solution

display('Done')
