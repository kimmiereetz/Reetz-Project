clear all
clc
% Setting up grid number of nodes N determines the fine-ness of the mesh
% the larger N the finer the mesh

N=60;   
% Based on the number of nodes the increments for dx dy and dt are
% determined
dx=2*pi/(N-1);
dt=dx^2/2;
% Based on the increment size for the time and space dimensions and their
% intervals the vectors for x, y, and time are created
x=0:dx:2*pi;
y=x;
t=0:dt:30;
% For ease of operating loops the number of elements in each x y and time
% vector are calculated
nx=length(x);
ts=length(t);
%%
% The constants determined from discritization are calculated befoe the
% loop to optimize code
a=dt/(2*dx^2);
b=(1+2*a);
c=(1-2*a);
% This mesh is required for plotting the surface solution 
[X,Y]=meshgrid(x,y);
%%
% The Dirichlet boundary conditions were reduced and are evaluated to 
% produce vectors.
%The Dirichlet boundary condition at x=0.
PHIab=cos(pi*y).*cosh(2*pi-y);
%The Dirichlet boundary condition at x=2?.
PSIab=y.^2.*sin(y/4);
%The Dirichlet boundary conditions are imposed on the grid boundaries while
% the initial condition is imposed on the interior points of the grid. 
U=[PHIab;zeros(nx-2,nx);PSIab];
Un=U;
%%
% For solving the implicit parts of the ADI method a tridiagonal equation
% matrix solver function was used
% Tridiag: Tridiagonal equation solver banded system
%   x = Tridiag(e,f,g,r): Tridiagonal system solver. 
% input: 
%   e = subdiagonal vector 
%   f = diagonal vector 
%   g = superdiagonal vector 
%   r = right hand side vector 
% output: 
%   x = solution vector

% For the first half step in time, the solution is implicit in x and the
% following subdiagonal, diagonal, and superdiagonal were determined for
% the corresponding tridiagonal matrix. Note the Dirichlet boundary
% conditions are imposed in this tridiagagonal system
% The right hand side vector is determined by the explict solution in y but
% the dimensions of the right hand side vector are preallocated preceding 
% this solution for code optimization

e=[-a*ones(1,N-2) 0];
f=[ 1 b*ones(1, N-2) 1];
g=[ 0 -a*ones(1, N-2)];
r=zeros(1,N);

% For the second half step in time, the solution is implicit in y and the
% following subdiagonal, diagonal, and superdiagonal were determined for
% the corresponding tridiagonal matrix. Note the Nuemann boundary
% conditions are imposed in this tridiagagonal system. Use of the ghost 
% nodes i.e. centered difference approximation is used for these boundary 
% conditions 
% The right hand side vector is determined by the explict solution in x but
% the dimensions of the right hand side vector are preallocated preceding 
% this solution for code optimization

e2=[-a*ones(1,N-2) -2*a];
f2=b*ones(N);
g2=[-2*a -a*ones(1,N-2)];
r2=zeros(1,N);
% The loops for solving the problem are ready to run. Time is the outer
% most loop 
% For the fisrt half step in time the columns of u(n+1/2) are transversed
% and determined using the explicit solution in y as the right hand side
% vector
for h=1:ts;
    for j=1:nx;
        for i=2:nx-1;
% r(i) is the explicit solution in y for the first half time step
% Nuemann boundary conditions are imposed in the explicit solution of y
% by use of the ghost nodes U(i,0) and U(i,N+1)
% where U(i,0)=U(i,2) and U(i,N+1)=U(i,N-1) i.e. the centered difference 
% approximation is used at the boundaries of y. Thus the discretization 
% reduces to r(i)=c*U(i,1)+2*a*U(i,2) at y=0 and 
% r(i)=2*a*U(i,nx-1)+c*U(i,nx) at y= 2pi          
            if j==1
                r(i)=c*U(i,1)+2*a*U(i,2);
            elseif j==nx
                r(i)=2*a*U(i,nx-1)+c*U(i,nx);
% Everywhere else the originally derived explicit side of the
            else
                r(i)=a*U(i,j-1)+c*U(i,j)+a*U(i,j+1);
            end
        end
% The Dirichlet boundary conditions are imposed on the explicit solution in
% y
        r=[U(1,j) r(2:nx-1) U(nx,j)];
% The implicit solution in x is determined using the tridiagonal equations
% solver in which the subdiagonal, diagonal, and superdiagonal where
% determined form discretization and the dirichlet boundary conditions for
% x
        x=tridiag(e,f,g,r);
% Finally u(n+1/2) is determined by the combination of the explicit
% solution in y and the implicit solution in x
        U(:,j)=x;
    end
    


%%
% For the second half step in time the rows of u(n+1) are transversed
% and determined using the explicit solution in x as the right hand side
% vector
    for i=2:nx-1;
        for j=1:nx;
% r(i) is the explicit solution in x for the second half time step            
             r2(j)=a*U(i-1,j)+c*U(i,j)+a*U(i+1,j);
        end
% The implicit solution in y is determined using the tridiagonal equations
% solver in which the subdiagonal, diagonal, and superdiagonal where
% determined form discretization and the Neumann boundary conditions for
% y        
        x2=tridiag(e2,f2,g2,r2);
% Finally u(n+1) is determined by the combination of the explicit
% solution in x and the implicit solution in x        
        U(i,:)=x2;
    end
% The centrally located grid node is tracked for use of steady state
% determination and diffusion confirmation

      center(h)=U(N/2,N/2);
% Throughout the time steps plots are generated to confirm diffusion
         figure (1)
     if mod(h,5)==0;
         surf(X,Y,U)
         drawnow
     end
        
end
% The centrally located grid node plotted versus time to confirm steady
% state and verify nature of diffusion
 figure (2)
  plot(t,center)
% The L2 Norm of the final U matrix is used for grid convergence. The L2
% Norms at Nodes of N 2N 4N 16N are calculted until the solution truncates 
% The minimum difference in L2 Norms is calculated to determine the sweet
% spot for the numerical solution
gi=norm(U,2);
