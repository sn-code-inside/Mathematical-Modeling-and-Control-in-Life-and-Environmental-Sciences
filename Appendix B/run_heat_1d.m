%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves the one-dimensional heat equation
% with homogeneous Neumann boundary conditions 
%
% u_t(x,t)-u_xx(x,t) = f(x,t)      in   (x_a,x_b) x (0,T)
% u_x(x_a,t) = u_x(x_b,t) = 0       
% u(x,0) = u_0(x)
%
% with P1 finite elements in space and 
% the Backward-Euler method in time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a = 0;
x_b = 1;
nx = input('nx = ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tend = 1;
dt = input('dt = ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIGHT HAND SIDE FUNCTION f(x,t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(x,t) zeros(size(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL DATUM u_0(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_0 = @(x) cos(pi*x.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXACT SOLUTION (uncomment the following lines for convergence tests) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uex = @(x,t) exp(-3*t).*cos(6*pi*x).^2;
f = @(x,t) -3*exp(-3*t).*cos(6*pi*x).^2+72*pi^2*exp(-3*t).*(cos(6*pi*x).^2-sin(6*pi*x).^2);
u_0 = @(x) uex(x,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MESH PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1 = nx+1;
h = (x_b-x_a)/nx;
nno = nx1;
nel = nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = [0:dt:Tend];
nt = length(tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NODE AND ELEM MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[node,elem] = get_mesh_P1_1d(x_a,x_b,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE STIFFNESS AND LUMPED MASS MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA,dmm] = get_mat_P1_1d(node,elem,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE THE ITERATION MATRIX MAT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mass = spdiags(dmm,0,nno,nno);
mat = 1/dt*Mass+AA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE CHOLESKY FACTORIZATION OF THE ITERATION MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR = chol(mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE THE SOLUTION VECTOR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu = zeros(nno,1);
uu = u_0(node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE THE OUTPUT SOLUTION MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_uu = zeros(nno,nt);
mat_uu(:,1) = uu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME LOOP START 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nt-1

   % extract current time step  
   t = tt(i+1);
 
   % assemble the RHS
   bb = f(node,t);   
   bb = 1/dt*uu+bb;
   bb = dmm.*bb;

   % solve the linear system
   yy = RR'\bb;
   uu = RR\yy;

   % update
   mat_uu(:,i+1) = uu;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nno
   for j=1:nt
      mat_x(i,j) = node(i);
      mat_t(i,j) = tt(j);
      mat_uex(i,j) = uex(node(i),tt(j));
   end
end

figure(1)
surf(mat_x',mat_t',mat_uu')
axis([0 1 0 1 0 1])
set(gca,'XTick',[0:0.2:1],'YTick',[0:0.2:1],'ZTick',[0:0.2:1],'Fontsize',30)
xlabel('x','Fontsize',30)
ylabel('t','Fontsize',30)
title('u_h(x,t)','Fontsize',30)

figure(2)
surf(mat_x',mat_t',mat_uex')
axis([0 1 0 1 0 1])
set(gca,'XTick',[0:0.2:1],'YTick',[0:0.2:1],'ZTick',[0:0.2:1],'Fontsize',30)
xlabel('x','Fontsize',30)
ylabel('t','Fontsize',30)
title('u(x,t)','Fontsize',30)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE ERROR (uncomment for convergence tests) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nt
   uu_ex = uex(node,tt(i));
   err(i) = norm(mat_uu(:,i)-uu_ex)/norm(uu_ex);
end
err = max(err);
disp(['error = ',num2str(err)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [node,elem] = get_mesh_P1_1d(x_a,x_b,nx)
% 
% This function constructs the data structure of a 1D mesh 
% for P1 finite elements
%
% INPUT
% 
% x_a := left vertex of the 1D interval domain
% x_b := right vertex of the 1D interval domain
% nx  := number of 1D finite elements (sub-intervals)
%
% OUTPUT
%
% node := node matrix ((nx+1)*1 matrix), 
%         containing the x-coordinates of the discretizations nodes 
% elem := elem (connectivity) matrix (nx*2), 
%         containing for each sub-interval, 
%         the indices of the two vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node,elem]=get_mesh_P1_1d(x_a,x_b,nx)

nx1 = nx+1;
h = (x_b-x_a)/nx;
nno = nx1;
nel = nx;

node = zeros(nno,1);
elem = zeros(nel,2);

for i = 1:nx1
   node(i,1) = x_a+(i-1)*h;
end

for i=1:nx
   elem(i,1) = i;
   elem(i,2) = i+1;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA,cc] = get_mat_P1_1d(node,elem,h)
%
% This function contsructs the stiffness and 
% mass matrices for 1D P1 finite elements
%
% INPUT
%
% node := node matrix ((nx+1)*1 matrix), 
%         containing the x-coordinates of the discretizations nodes
% elem := elem (connectivity) matrix (nx*2), 
%         containing for each sub-interval, 
%         the indices of the two vertices
% h    := mesh size (length of the generic sub-interval) 
%
% OUTPUT
%
% AA := stiffness matrix
% cc := diagonal mass matrix, obtained using mass lumping 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AA,cc] = get_mat_P1_1d(node,elem,h) 

% get number of nodes (nno) and elements (nel)
nno=size(node,1);
nel=size(elem,1);

% initiliaze AA and cc
AA = sparse(nno,nno);
cc = zeros(nno,1);

% compute local stiffness and mass matrices
[AA_loc,cc_loc] = get_mat_P1_1d_loc(h);

% assembling loop on the elements
for i=1:nel

   % extract the indices of the vertices of i-th element
   kpm(1) = elem(i,1);
   kpm(2) = elem(i,2);

   % add local contributions
   AA(kpm,kpm) = AA(kpm,kpm)+AA_loc;
   cc(kpm) = cc(kpm)+cc_loc;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA_loc,cc_loc] = get_mat_P1_1d_loc(h)
%
% This function constructs the local stiffness and 
% mass matrices for 1D P1 finite elements
%
% INPUT
%
% h    := mesh size (length of the generic sub-interval)
%
% OUTPUT
%
% AA_loc := local stiffness matrix
% cc_loc := local diagonal mass matrix, obtained using mass lumping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AA_loc,cc_loc] = get_mat_P1_1d_loc(h) 

AA_loc = zeros(2,2);
cc_loc = zeros(2,1);

AA_loc(1,1) = 1/h;
AA_loc(2,2) = 1/h;
AA_loc(1,2) = -1/h;
AA_loc(2,1) = -1/h;

cc_loc(1) = h/2;
cc_loc(2) = h/2;

end
