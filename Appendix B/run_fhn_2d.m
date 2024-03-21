%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves the FitzHugh-Nagumo system in 2D 
%
% u_t - sigma laplacian(u) + f(u,w) = I_app  in \Omega=(0,1)^2 x (0,T)
% w_t + g(u,w) = 0                           in \Omega=(0,1)^2 x (0,T)
% sigma grad(u) n = 0                        on \partial(\Omega) x (0,T) 
% u(0) = u_0                                 in \Omega=(0,1)^2
% w(0) = w_0                                 in \Omega=(0,1)^2
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
nx = 100;
ny = nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tend = input('Tend = ');
dt = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 2e-3;
b = 5;
beta = 0.1;
delta = 1;
c = 1;
gamma = 0.25;
e = 0.1;
iap = 5;
rad = 0.04;
T_sti = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS f(u,w) and g(u,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(u,w) -b*u.*(u-beta).*(delta-u)+c*w;
g = @(u,w) -e*(u-gamma*w);

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
% SET TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = [0:dt:Tend];
nt = length(tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NODE AND ELEM MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[node,elem] = get_mesh_Q1_2d(x_a,x_b,y_a,y_b,nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE STIFFNESS AND LUMPED MASS MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA,dmm] = get_mat_Q1_2d(node,elem,h);
AA = sigma*AA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE THE ITERATION MATRIX MAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mass = spdiags(dmm,0,nno,nno);
mat = [1/dt*Mass+AA      sparse(nno,nno);
       sparse(nno,nno)   1/dt*Mass      ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE CHOLESKY FACTORIZATION OF THE ITERATION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR = chol(mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE SOLUTION VECTORS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu = zeros(nno,1);
ww = zeros(nno,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE APPLIED CURRENT VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iapp = zeros(nno,1);
for i = 1:nno
   x = node(i,1);
   y = node(i,2);
   r = sqrt((x-0.5)^2+(y-0.5)^2);
   if(r<=rad);
     Iapp(i) = iap;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE OUTPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iind = ny/4*nx1+nx/4+1;
uu_1(1) = 0;
ww_1(1) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME LOOP START
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nt-1

   % extract current time step
   t = tt(i+1);

   % assemble the RHS
   bb_u = -f(uu,ww);
   if(t<=T_sti)
    bb_u = bb_u+Iapp;
   end
   bb_u = 1/dt*uu+bb_u;
   bb_u = dmm.*bb_u;
   bb_w = -g(uu,ww);
   bb_w = 1/dt*ww+bb_w;
   bb_w = dmm.*bb_w;
   bb = [bb_u;bb_w];

   % solve the linear system
   yy = RR'\bb;
   uw = RR\yy;

   % update
   uu = uw(1:nno);
   ww = uw(nno+1:end);

   % save u and w in a sample point
   uu_1(i+1) = uu(iind);
   ww_1(i+1) = ww(iind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ny1
   for i=1:nx1
      iind = (j-1)*nx1+i;
      mat_x(j,i) = node(iind,1);
      mat_y(j,i) = node(iind,2);
      mat_uu(j,i) = uu(iind,1);
   end
end
han = subplot(1,1,1);
surf(mat_x,mat_y,mat_uu)
axis([0 1 0 1 -0.5 1])
set(han,'CLim',[-0.5 1]);
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1],'ZTick',[-0.5:0.5:1],'Fontsize',50)
title(['t = ',num2str(Tend)],'Fontsize',50)

figure
subplot(1,2,1)
plot(tt,uu_1,'k','Linewidth',2)
axis([0 Tend -0.5 1])
set(gca,'XTick',[0:5:Tend],'YTick',[-0.5:0.5:1],'Fontsize',30)
title('u(x,t)','Fontsize',30)
subplot(1,2,2)
plot(tt,ww_1,'k','Linewidth',2)
axis([0 Tend 0 0.8])
set(gca,'XTick',[0:5:Tend],'YTick',[0:0.2:0.8],'Fontsize',30)
title('w(x,t)','Fontsize',30)


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
