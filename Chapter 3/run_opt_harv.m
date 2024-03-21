%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM RUN_OPT_HARV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program computes the regional control in the optimal harvesting problem,
% equations (3.13), (3.14), (3.18),
% according to the conceptual algorithm reported in section 3.5.1,
% on the square domain [x_a,x_b] x [y_a,y_b]
% and the time interval [0,Tend].
% With the current set of parameters, it reproduces in particular
% the results reported in the bottom right panel of Figure 3.1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1;
L = 1;
epsilon = 1;
c1 = 0.6;
c2 = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a = 0;
x_b = 1;
y_a = 0;
y_b = 1;
nx = 64;         % number of mesh elements in x-direction
ny = 64;         % number of mesh elements in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tend = 1;
dt = 0.05;       % time step size
dth = 0.00001;   % fictitious time step size for the solution of equation (1.19)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOLERANCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol_1 = 1e-3;    % this is epsilon_1 in the conceptual algorithm 
tol_2 = 1e-3;    % this is epsilon_2 in the conceptual algorithm 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MODEL PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param(1) = d;
param(2) = L;
param(3) = epsilon;
param(4) = c1;
param(5) = c2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MESH PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1 = nx+1;
ny1 = ny+1;
h = (x_b-x_a)/nx;    % mesh size
nno = nx1*ny1;       % number of nodes
nel = nx*ny;         % number of elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET TIME PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = [0:dt:Tend];    % time steps
nt = length(tt);     % number of time steps 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NODE AND ELEM MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[node,elem] = get_mesh_Q1_2d(x_a,x_b,y_a,y_b,nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE STIFFNESS AND LUMPED MASS MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA,dmm] = get_mat_Q1_2d(node,elem,h);

AA_1 = d*AA;
AA_2 = d*AA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FINITE ELEMENT STRUCTURE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fem = struct('node',node,'elem',elem,'AA_1',AA_1,'AA_2',AA_2,'dmm',dmm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIAL CONDITIONS FOR STATE VARIABLES r(x,t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rr0 = zeros(nno,1);
for i=1:nno
   x = node(i,1);
   y = node(i,2);

   rr0(i) = 1/(2*pi)*exp(-(x^2+y^2)/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE CONTROLS phi(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = zeros(nno,1);
for i=1:nno
   x = node(i,1);
   y = node(i,2);

   phi(i) = (-x+0.25)*(-y+0.5)+1/4*sin(x-0.4)*sin(0.25-y);
%   phi(i) = sin(5*pi*x)*sin(5*pi*y);
end
phi0 = phi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OPTIMIZATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = 1e+6;
op_it = 0;

while(op_it<1000) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE BACKWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('solve the backward problem')
  mpp = solve_pback(phi,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE COST FUNCTIONAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('compute the cost functional')
  J_ = comp_costJ(rr0,mpp(:,1),phi,param,fem);
  disp(['J_ = ',num2str(J_)]);

  err_J = abs(J_-J);
  
  if(err_J<tol_1 || J_>=J)
    break
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE FORWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('solve the forward problem')
  mrr = solve_pforw(rr0,mpp,phi,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE THE CONTROL PHI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('update the control phi')
  phi_ = update_phi(phi,mpp,mrr,param,fem,tt,dth);
  err_phi = norm(phi_-phi);
  disp(['err phi = ',num2str(err_phi)]);

  if(err_phi<tol_2)
    break
  end

  phi = phi_;
  J = J_;
  op_it =op_it+1;
  disp(['op_it = ',num2str(op_it)]);
  Jv(op_it) = J;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ny1
   for i=1:nx1
      iind = (j-1)*nx1+i;
      mat_pp(j,i) = mpp(iind,1);
      mat_rr(j,i) = mrr(iind,end);
      mat_phi0(j,i) = phi0(iind);
      mat_phi(j,i) = phi(iind);
   end
end

figure(1)
surf(mat_pp)
title('p(x,y)','Fontsize',30)

figure(2)
surf(mat_rr)
title('r(x,y)','Fontsize',30)

figure(3)
contourf(mat_phi0,[0 0])
axis equal
axis tight

figure(4)
contourf(mat_phi,[0 0])
axis equal
axis tight
set(gca,'XTick',[],'YTick',[])

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function mrr = solve_pforw(rr0,mpp,phi,param,fem,dt,tt)
%
% This functions solves the forward problem, i.e. the equation for r(x,t) 
%
% INPUT
%
% rr0   := initial condition for the state variable r(x,t) 
% mpp   := spatio-temporal evolution of p(x,t) 
% mphi  := distribution of control phi(x) 
% param := parameters
% fem   := finite element structure
% dt    := timestep size
% tt    := vector of timesteps
%
% OUTPUT
%
% mrr   := spatio-temporal evolution of the state variable r(x,t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mrr = solve_pforw(rr0,mpp,phi,param,fem,dt,tt) 

% get the number of nodes 
nno = size(rr0,1);

% set time loop parameters
nt = length(tt);

% set initial conditions
mrr = zeros(nno,nt);
mrr(:,1) = rr0;

% assemble iteration matrix
Mass = spdiags(fem.dmm,0,nno,nno);
mat = 1/dt*Mass+fem.AA_2;

% time loop start
for k=1:nt-1
   t = tt(k+1);

   % compute and integrate reaction term of equation for k(x,t)
   rt_1 = rt1_fun(mrr(:,k),mpp(:,k),phi,param,fem);
   bb_1 = fem.dmm.*(1/dt*mrr(:,k)+rt_1);

   % solve the linear system and update
   mrr(:,k+1)=mat\bb_1;
end

clear rt_1;
clear Mass;
clear mat;
clear bb_1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rt_1 = rt1_fun(rr,pp,phi,param,fem); 
%
% This functions constructs the reaction term of the equation for r(x,t) 
%
% INPUT
%
% rr    := current state variable r(x,t) distribution
% pp    := current state variable p(x,t) distribution
% phi   := control phi(x) distribution
% param := parameters
% fem   := finite element structure
%
% OUTPUT
%
% rt_1  := reaction term of the equation for r(x,t) evaluated in each mesh node 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt_1 = rt1_fun(rr,pp,phi,param,fem); 

% get the number of nodes
nno = size(rr,1);

% initialize the output
rt_1 = zeros(nno,1);

% extract useful parameters
L = param(2);
epsilon = param(3);

% compute eta and H
x = fem.node(:,1);
y = fem.node(:,2);
eta = eta_fun(x,y);
Hphi = H_fun(phi,epsilon);
Hpp = H_fun(1+pp,epsilon);
dpp = delta_fun(1+pp,epsilon);

% compute the reaction term
rt_1 = (eta-L*Hphi.*(Hpp+(1+pp).*dpp)).*rr;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function mpp = solve_pback(phi,param,fem,dt,tt)
%
% This functions solves the backward problem, i.e. the equation for p(x,t) 
%
% INPUT
%
% phi   := control phi(x)
% param := parameters
% fem   := finite element structure
% dt    := timestep size
% tt    := vector of timesteps
%
% OUTPUT
%
% mpp   := spatio-temporal evolution of p(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpp = solve_pback(phi,param,fem,dt,tt)

% get the number of nodes
nno = size(phi,1);

% set time loop parameters
nt = length(tt);

% set initial conditions
mpp = zeros(nno,nt);

% assemble iteration matrix
Mass = spdiags(fem.dmm,0,nno,nno);
mat = 1/dt*Mass+fem.AA_1;

% time loop start
for k=1:nt-1
   t = tt(end-k);

   % compute and integrate reaction term of equation for lambda_k(x,t)
   rt_1 = rtb1_fun(mpp(:,nt-k+1),phi,param,fem);
   bb_1 = fem.dmm.*(1/dt*mpp(:,nt-k+1)-rt_1);

   % solve the linear system and update
   mpp(:,nt-k)=mat\bb_1;
end

clear rt_1;
clear Mass;
clear mat;
clear bb_1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rt_1 = rtb1_fun(pp,phi,param,fem)
%
% This functions constructs the reaction term of the equation for p(x,t) 
%
% INPUT
%
% pp    := current variable p(x,t) distribution
% phi   := control phi(x) distribution
% param := parameters
% fem   := finite element structure
%
% OUTPUT
%
% rt_1  := reaction term of the first equation evaluated in each mesh node 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt_1 = rtb1_fun(pp,phi,param,fem)

% get the number of nodes
nno = size(pp,1);

% initialize the output
rt_1 = zeros(nno,1);

% extract useful parameters
L = param(2);
epsilon = param(3);

% compute eta and H
x = fem.node(:,1);
y = fem.node(:,2);
eta = eta_fun(x,y);
Hphi = H_fun(phi,epsilon); 
Hpp = H_fun(1+pp,epsilon);

% compute the reaction term
rt_1 = -eta.*pp+L*Hphi.*(1+pp).*Hpp;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Hu = H_fun(u,epsilon)
%
% This functions evaluates H_\epsilon(u) 
%
% INPUT
%
% u         := current u distribution
% epsilon   := epsilon parameter 
%
% OUTPUT
%
% Hu        := H_\epsilon(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hu = H_fun(u,epsilon)

Hu = 1/2+(1+2/pi*atan(u/epsilon));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function eta = eta_fun(x,y)
%
% This functions evaluates eta(x,y)
%
% INPUT
%
% x   := x-coordinate 
% y   := y-coordinate 
%
% OUTPUT
%
% eta := eta(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eta = eta_fun(x,y)

eta = 3*ones(size(x));
%eta = exp(1.5*(x+y)); 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function du = delta_fun(u,epsilon)
%
% This functions evaluates \delta_\epsilon(u)
%
% INPUT
%
% u         := current u distribution
% epsilon   := epsilon parameter
%
% OUTPUT
%
% du        := \delta_\epsilon(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function du = delta_fun(u,epsilon)

du = epsilon./(pi.*(epsilon^2+u.^2));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ddu = ddelta_fun(u,epsilon)
% 
% This functions evaluates \delta^\prime_\epsilon(u)
%
% INPUT
%
% u         := current u distribution
% epsilon   := epsilon parameter
%
% OUTPUT
% 
% ddu        := \delta^\prime_\epsilon(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddu = ddelta_fun(u,epsilon)

ddu = -2*epsilon/pi*u./(epsilon^2+u.^2).^2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Aphi = comp_Aphi(node,elem,phi,param)
%
% This function computes the non-linear diffusion term of equation (1.19) 
%
% INPUT
%
% node   := node matrix, containing the x- and y-coordinates of 
%           the discretizations nodes
% elem   := elem (connectivity) matrix, containing for each 
%           element, the indices of the four vertices
% phi    := control phi distribution 
% param := parameters
%
% OUTPUT
%
% Aphi   := the non-linear diffusion term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aphi = comp_Aphi(node,elem,phi,param) 

% get the number of nodes and elements
nno = size(node,1);
nel = size(elem,1);

% initialize the output
Aphi = zeros(nno,1);

% extract useful parameters
epsilon = param(3);

% get the Q1 basis functions and gradients
pphi = get_phi();
[pxs,pys] = get_dphi();

% loop on the elements
for i = 1:nel
   l1 = elem(i,1);
   l2 = elem(i,2);
   l3 = elem(i,3);
   l4 = elem(i,4);

   kpm(1) = l1;
   kpm(2) = l2;
   kpm(3) = l3;
   kpm(4) = l4;	 

   for k=1:4
      for j=1:2
         c(j,k) = node(kpm(k),j);
      end
   end

   phi_loc = phi(kpm);

   Aphi_loc = comp_Aphi_loc(c,pphi,pxs,pys,phi_loc,epsilon);
	
   Aphi(kpm) = Aphi(kpm)+Aphi_loc;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Aphi_loc = comp_Aphi_loc(c,pphi,pxs,pys,phi_loc,epsilon)
%
% This function computes the local (to the mesh element) 
% non-linear diffusion term of equation (1.19) 
%
% INPUT
%
% c         := coordinates of the vertices of the element
% pphi      := Q1 basis functions
% pxs,pys   := gradients of Q1 basis functions
% phi_loc   := components of control phi local to the element 
% epsilon   := epsilon parameter
%
% OUTPUT
%
% Aphi_loc  := the local non-linear diffusion term 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aphi_loc = comp_Aphi_loc(c,pphi,pxs,pys,phi_loc,epsilon)

for l = 1:4
   xx(l) = c(1,l);
   yy(l) = c(2,l);
end
Aphi_loc = zeros(4,1);

for l = 1:4
   for i = 1:2
      for j = 1:2
         t(i,j) = 0;
         s(i,j) = 0;
      end
   end
   for k = 1:4
      t(1,1) = t(1,1)+xx(k)*pxs(k,l);
      t(1,2) = t(1,2)+xx(k)*pys(k,l);
      t(2,1) = t(2,1)+yy(k)*pxs(k,l);
      t(2,2) = t(2,2)+yy(k)*pys(k,l);
   end
   [d,s] = det_inv(t);
   des(l) = abs(d);
   for k=1:4
      px(k,l) = s(1,1)*pxs(k,l)+s(1,2)*pys(k,l);
      py(k,l) = s(2,1)*pxs(k,l)+s(2,2)*pys(k,l);
   end
end

for l=1:4
   dphi(1) = dot(phi_loc,px(:,l));
   dphi(2) = dot(phi_loc,py(:,l));
   ndphi = norm(dphi);
   if(ndphi<1e-8)
     ndphi = 1e-8;
   end
   del = delta_fun(phi_loc(l),epsilon);
   ddel = ddelta_fun(phi_loc(l),epsilon);
   for i=1:4
      phi(1) = px(i,l);
      phi(2) = py(i,l);
      pphil = pphi(i,l);
      ss = ddel*dot(dphi,dphi)/ndphi*pphil+del*dot(dphi,phi)/ndphi;
      Aphi_loc(i) = Aphi_loc(i)+ss*des(l)/4;
   end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function phi = get_phi()
%
% This function computes the Q1 basis functions at the quadrature points 
%
% OUTPUT
%
% phi   := Q1 basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = get_phi()

phi(1,1) = 1;
phi(1,2) = 0;
phi(1,3) = 0;
phi(1,4) = 0;
 
phi(2,1) = 0;
phi(2,2) = 1;
phi(2,3) = 0;
phi(2,4) = 0;

phi(3,1) = 0;
phi(3,2) = 0;
phi(3,3) = 1;
phi(3,4) = 0;

phi(4,1) = 0;
phi(4,2) = 0;
phi(4,3) = 0;
phi(4,4) = 1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dphidx,dphidy] = get_dphi()
%
% This function computes the gradients of the Q1 basis 
% functions at the quadrature points 
%
% OUTPUT
%
% dphidx,dphidy   := gradients of Q1 basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dphidx,dphidy] = get_dphi()

dphidx(1,1) = -1;
dphidx(1,2) = -1;
dphidx(1,3) = 0;
dphidx(1,4) = 0;

dphidx(2,1) = 1;
dphidx(2,2) = 1;
dphidx(2,3) = 0;
dphidx(2,4) = 0;

dphidx(3,1) = 0;
dphidx(3,2) = 0;
dphidx(3,3) = -1;
dphidx(3,4) = -1;

dphidx(4,1) = 0;
dphidx(4,2) = 0;
dphidx(4,3) = 1;
dphidx(4,4) = 1;

dphidy(1,1) = -1;
dphidy(1,2) = 0;
dphidy(1,3) = -1;
dphidy(1,4) = 0;

dphidy(2,1) = 0;
dphidy(2,2) = -1;
dphidy(2,3) = 0;
dphidy(2,4) = -1;

dphidy(3,1) = 1;
dphidy(3,2) = 0;
dphidy(3,3) = 1;
dphidy(3,4) = 0;

dphidy(4,1) = 0;
dphidy(4,2) = 1;
dphidy(4,3) = 0;
dphidy(4,4) = 1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [d,S] = det_inv(t)
%
% This function computes the determinant and 
% the inverse of a square matrix 
%
% INPUT
%
% t    := square matrix 
%
% OUTPUT
%
% d    := determinant of t
% S    := inverse of t 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,S] = det_inv(t)

d = det(t);
S = inv(t)';

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function phi_ = update_phi(phi,mpp,mrr,param,fem,tt)
%
% This function updates the contro phi by solving equation (1.19) 
%
% INPUT
%
% phi   := distribution of control phi(x)
% mpp   := spatio-temporal evolution of p(x,t)
% mrr   := spatio-temporal evolution of r(x,t)
% param := parameters
% fem   := finite element structure
% tt    := vector of timesteps
%
% OUTPUT
%
% phi_  := updated distribution of control phi(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi_ = update_phi(phi,mpp,mrr,param,fem,tt,dth) 

% get the number of nodes
nno=size(phi,1);

% get the number of time steps
nt = length(tt);

% extract useful parameters
L = param(2);
epsilon = param(3);
c1 = param(4);
c2 = param(5);

% compute Ft_1, the reaction term of equation (1.19)
Ft_1 = zeros(nno,1);

for j=1:nno
   iF = zeros(nt,1);
   for i = 1:nt
      Hphi_l = H_fun(1+mpp(j,i),epsilon);
      iF(i) = (1+mpp(j,i))*Hphi_l*mrr(j,i); 
   end
   Ft_1(j)=trapz(tt,iF);
end

% start theta loop
for k = 1:100
   % compute \delta_\epsilon(phi)
   dphi = delta_fun(phi,epsilon);

   % compute Aphi, the diffusion term of equation 1.19 
   Aphi = comp_Aphi(fem.node,fem.elem,phi,param);

   % update phi
   bb = -c1*Aphi+fem.dmm.*dphi.*(-c2+L*Ft_1);
   phi_ = phi+dth*bb./fem.dmm;
   
   phi = phi_;
end

clear bb;
clear Ft_1;
clear Aphi;
clear dphi;
clear iF;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function J = comp_costJ(rr0,pp0,phi,param,fem)
%
% This function computes the cost functional 
%
% INPUT
%
% rr0   := distribution of r(x,t) at t = 0
% pp0   := distribution of p(x,t) at y = 0
% phi   := distribution of control phi(x)
% param := parameters
% fem   := finite element structure
%
% OUTPUT
%
% J     := cost functional 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = comp_costJ(rr0,pp0,phi,param,fem) 

% get the number of nodes
nno=size(phi,1);

% extract useful parameters
epsilon = param(3);
c1 = param(4);
c2 = param(5);

% compute J_1
J_1 = dot(fem.dmm,rr0.*pp0);

% compute J_2
J_2 = comp_J_2(fem.node,fem.elem,phi,epsilon);
J_2 = c1*J_2;

% compute J_3
Hphi = H_fun(phi,epsilon);
J_3 = c2*dot(fem.dmm,Hphi);

% compute the total cost J
J = J_1+J_2+J_3;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function J_2 = comp_J_2(node,elem,phi,epsilon)
%
% This function computes the second term of the cost functional 
%
% INPUT
%
% node      := node matrix, containing the x- and y-coordinates 
%              of the discretizations nodes
% elem      := elem (connectivity) matrix, containing for each 
%              element, the indices of the four vertices
% phi       := control phi distribution
% epsilon   := epsilon parameter
%
% OUTPUT
%
% J_2       := second term of the cost functional 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J_2 = comp_J_2(node,elem,phi,epsilon) 

% get the number of nodes and elements
nno = size(node,1);
nel = size(elem,1);

% initialize the output
J_2 = 0;

% get the Q1 basis functions gradients
[pxs,pys] = get_dphi();

% loop on the elements
for i = 1:nel
   l1 = elem(i,1);
   l2 = elem(i,2);
   l3 = elem(i,3);
   l4 = elem(i,4);

   kpm(1) = l1;
   kpm(2) = l2;
   kpm(3) = l3;
   kpm(4) = l4;

   for k=1:4
      for j=1:2
         c(j,k) = node(kpm(k),j);
      end
   end

   phi_loc = phi(kpm);

   J_2_loc = comp_J_2_loc(c,pxs,pys,phi_loc,epsilon);

   J_2 = J_2+J_2_loc;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function J_2_loc = comp_J_2_loc(c,pxs,pys,phi_loc,epsilon)
%
% This function computes the local part of the second term of 
% the cost functional 
%
% INPUT
%
% c         := coordinates of the vertices of the element
% pxs,pys   := gradients of Q1 basis functions
% phi_loc   := components of control phi local to the element
% epsilon   := epsilon parameter
%
% OUTPUT
%
% J_2_loc   := local part of the second term of the cost functional 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J_2_loc = comp_J_2_loc(c,pxs,pys,phi_loc,epsilon)

for l = 1:4
   xx(l) = c(1,l);
   yy(l) = c(2,l);
end
J_2_loc = 0;

for l = 1:4
   for i = 1:2
      for j = 1:2
         t(i,j) = 0;
         s(i,j) = 0;
      end
   end
   for k = 1:4
      t(1,1) = t(1,1)+xx(k)*pxs(k,l);
      t(1,2) = t(1,2)+xx(k)*pys(k,l);
      t(2,1) = t(2,1)+yy(k)*pxs(k,l);
      t(2,2) = t(2,2)+yy(k)*pys(k,l);
   end
   [d,s] = det_inv(t);
   des(l) = abs(d);
   for k=1:4
      px(k,l) = s(1,1)*pxs(k,l)+s(1,2)*pys(k,l);
      py(k,l) = s(2,1)*pxs(k,l)+s(2,2)*pys(k,l);
   end
end

for l=1:4
   dphi(1) = dot(phi_loc,px(:,l));
   dphi(2) = dot(phi_loc,py(:,l));
   ndphi = norm(dphi);
   del = delta_fun(phi_loc(l),epsilon);
   ss = del*ndphi;
   J_2_loc = J_2_loc+ss*des(l)/4;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% node := node matrix (((nx+1)*(ny+1))*2 matrix), 
%         containing the x- and y-coordinates of the discretizations nodes
% elem := elem (connectivity) matrix ((nx*ny)*4), containing 
%         for each element, the indices of the four vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA,cc] = get_mat_Q1_2d(node,elem,h)
%
% This functions constructs the stiffness and 
% mass matrices for 2D Q1 finite elements
%
% INPUT
%
% node := node matrix, containing the x- and y-coordinates of 
%         the discretizations nodes
% elem := elem (connectivity) matrix, containing for each 
%         element, the indices of the four vertices
% h    := mesh size (length of the generic sub-interval)
%
% OUTPUT
%
% AA := stiffness matrix
% cc := diagonal mass matrix, obtained using mass lumping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA_loc,cc_loc] = get_mat_Q1_2d_loc(h)
%
% This functions constructs the local stiffness and 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

