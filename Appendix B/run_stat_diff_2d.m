%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves the stationary diffusion problem
%
% -laplacian(u(x))  = f(x)      in \Omega=(0,1)^2
% u = g                         on \partial(\Omega) 
%
% with Q1 finite elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a = 0;
x_b = 1;
y_a = 0;
y_b = 1;
nx = input('nx = ');
ny = nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIGHT HAND SIDE FUNCTION f(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(x,y) -8*x+2*y.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXACT SOLUTION (uncomment the following lines for convergence tests) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uex = @(x,y) sin(8*pi*x).*y.^2-exp(x-4*y);
f = @(x,y) 64*pi^2*sin(8*pi*x).*y.^2+exp(x-4*y)-2*sin(8*pi*x)+16*exp(x-4*y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MESH PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1 = nx+1;
ny1 = ny+1;
h = (x_b-x_a)/nx;
nno = nx1*ny1;
nel = nx*nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NODE AND ELEM MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[node,elem] = get_mesh_Q1_2d(x_a,x_b,y_a,y_b,nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE STIFFNESS AND LUMPED MASS MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA,dmm] = get_mat_Q1_2d(node,elem,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET DIRICHLET AND INTERIOR DOFS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = 0;
ddofs = [];
for i = 1:nx1
   kk = kk+1;
   ddofs(kk) = i;
end
for i = 1:ny-1
   kk = kk+1;
   ddofs(kk) = i*nx1+1;
   kk = kk+1;
   ddofs(kk) = (i+1)*nx1;
end
for i = 1:nx1
   kk = kk+1;
   ddofs(kk) = ny*nx1+i;
end
idofs = setdiff([1:nno],ddofs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE SOLUTION VECTOR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu = zeros(nno,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET DIRICHLET BOUNDARY CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(ddofs)
   x = node(ddofs(i),1);
   y = node(ddofs(i),2);
   uu(ddofs(i)) = uex(x,y);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE THE RHS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bb = dmm.*f(node(:,1),node(:,2));
bb = bb-AA*uu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE LINEAR SYSTEM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu(idofs) = AA(idofs,idofs)\bb(idofs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE ERROR IN l2 NORM (uncomment for convergence tests) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu_ex = uex(node(:,1),node(:,2));
err = norm(uu-uu_ex)/norm(uu_ex);
disp(['l2 error = ',num2str(err)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ny1
   for i=1:nx1
      iind = (j-1)*nx1+i;
      mat_x(j,i) = node(iind,1);
      mat_y(j,i) = node(iind,2);
      mat_uu(j,i) = uu(iind);
      mat_uex(j,i) = uu_ex(iind);
   end
end

figure(1)
surf(mat_x,mat_y,mat_uu)
axis([0 1 0 1 -3 1])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1],'ZTick',[-3:1:1],'Fontsize',50)
title('u_h(x,y)','Fontsize',30)

figure(2)
surf(mat_x,mat_y,mat_uex)
axis([0 1 0 1 -3 1])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1],'ZTick',[-3:1:1],'Fontsize',50)
title('exact solution','Fontsize',30)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [node,elem] = get_mesh_Q1_2d(x_a,x_b,y_a,y_b,nx,ny)
%
% This functions constructs the data structure of a 2D mesh 
% for Q1 finite elements
%
% INPUT
%
% x_a := left x-coordinate of the 2D rectangular domain [x_a,x_b] x [y_a,y_b]
% x_b := right x-coordinate of the 2D rectangular domain [x_a,x_b] x [y_a,y_b] 
% y_a := left y-coordinate of the 2D rectangular domain [x_a,x_b] x [y_a,y_b]
% y_b := right y-coordinate of the 2D rectangular domain [x_a,x_b] x [y_a,y_b]
% nx  := number of elements (sub-intervals) along the x-direction
% ny  := number of elements (sub-intervals) along the y-direction
%
% OUTPUT
%
% node := node matrix (((nx+1)*(ny+1))*2 matrix), containing the x- and 
%         y-coordinates of the discretizations nodes
% elem := elem (connectivity) matrix ((nx*ny)*4), containing for each element, 
%         the indices of the four vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node,elem] = get_mesh_Q1_2d(x_a,x_b,y_a,y_b,nx,ny)

nx1 = nx+1;
ny1 = ny+1;
hx = (x_b-x_a)/nx;
hy = (y_b-y_a)/ny;
nno = nx1*ny1;
nel = nx*ny;

node = zeros(nno,2);
elem = zeros(nel,4);

for j=1:ny1
   for i=1:nx1
      node((j-1)*nx1+i,1)=x_a+(i-1)*hx;
      node((j-1)*nx1+i,2)=y_a+(j-1)*hy;
   end
end

for j=1:ny 
   for i=1:nx
      elem((j-1)*nx+i,1) = (j-1)*nx1+i;
      elem((j-1)*nx+i,2) = (j-1)*nx1+i+1;
      elem((j-1)*nx+i,3) = j*nx1+i;
      elem((j-1)*nx+i,4) = j*nx1+i+1;
   end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA,cc] = get_mat_Q1_2d(node,elem,h)
%
% This functions contsructs the stiffness and 
% mass matrices for 1D P1 finite elements
%
% INPUT
%
% node := node matrix, containing the x- and y-coordinates 
%         of the discretizations nodes
% elem := elem (connectivity) matrix, containing for each element, 
%         the indices of the four vertices
% h    := mesh size (length of the generic sub-interval)
%
% OUTPUT
%
% AA := stiffness matrix
% cc := diagonal mass matrix, obtained using mass lumping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AA,cc] = get_mat_Q1_2d(node,elem,h) 

nno = size(node,1);
nel = size(elem,1);

AA = sparse(nno,nno);
cc = zeros(nno,1);

[AA_loc,cc_loc] = get_mat_Q1_2d_loc(h);

for i=1:nel
   l1 = elem(i,1);
   l2 = elem(i,2);
   l3 = elem(i,3);
   l4 = elem(i,4);

   kpm(1) = l1;
   kpm(2) = l2;
   kpm(3) = l3;
   kpm(4) = l4;	 

   AA(kpm,kpm) = AA(kpm,kpm)+AA_loc;
   cc(kpm) = cc(kpm)+cc_loc;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA_loc,cc_loc] = get_mat_Q1_2d_loc(h)
%
% This functions contsructs the local stiffness and 
% mass matrices for 2D Q1 finite elements
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
function [AA_loc,cc_loc] = get_mat_Q1_2d_loc(h) 

AA_loc=zeros(4,4);
cc_loc=zeros(4,1);

AA_loc(1,1) = 2/3;
AA_loc(2,2) = 2/3;
AA_loc(3,3) = 2/3;
AA_loc(4,4) = 2/3;

AA_loc(1,2) = -1/6;
AA_loc(2,1) = -1/6;
AA_loc(3,4) = -1/6;
AA_loc(4,3) = -1/6;

AA_loc(1,3) = -1/6;
AA_loc(3,1) = -1/6;
AA_loc(2,4) = -1/6;
AA_loc(4,2) = -1/6;

AA_loc(1,4) = -1/3;
AA_loc(4,1) = -1/3;
AA_loc(2,3) = -1/3;
AA_loc(3,2) = -1/3;

for i=1:4
   cc_loc(i)=h^2/4;
end

end
