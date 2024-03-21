%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM RUN_XYLELLA_PRED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program solves the controlled system presented in Chapter 6,
% modeling the spread of Xylella epidemic and its control with predators,
% on the rectangular domain [x_a,x_b] x [y_a,y_b] 
% and the time interval [0,Tend].
% With the current set of parameters, it reproduces in particular
% the results of case B22.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1e-4;
r = 200;
chi = 0.01;
n = 0.98;
q = 0.5;
l = 0.01;
C = 100;
zeta = 0.2;
b = 0.05;
alpha = 0.3;
mu = 0.9;
p = 1;
c = 0.3;
r_2 = 2/3*r;
d_Z = 0.1*d;
delta_1 = 80;
delta_2 = 0.3;
%delta_1 = 1500;
%delta_2 = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a = 0;
x_b = 4;
y_a = 0;
y_b = 0.8;
nx = 200;
ny = 40;
h = (x_b-x_a)/nx;       % mesh size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tend = 20;
dt = 0.0001;            % time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_21 = 0;
gamma_1 = 0;
gamma_23 = 0;
gamma_22 = 0;
gamma_3 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS VECTOR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param(1) = d;
param(2) = r;
param(3) = chi;
param(4) = n;
param(5) = q;
param(6) = l;
param(7) = C;
param(8) = zeta;
param(9) = b;
param(10) = alpha;
param(11) = mu;
param(12) = p;
param(13) = c;
param(14) = r_2;
param(15) = d_Z;
param(16) = delta_1;
param(17) = delta_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET CONTROLS VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
control(1) = gamma_21;
control(2) = gamma_1;
control(3) = gamma_23;
control(4) = gamma_22;
control(5) = gamma_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1 = nx+1;
ny1 = ny+1;
nno = nx1*ny1;      % number of nodes
nel = nx*ny;        % number of Q1 finite elements
dofs = 5*nno;       % number of degrees of freedom 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = [0:dt:Tend];   % vector of timesteps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NODE AND ELEM MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[node,elem] = get_mesh_Q1_2d(x_a,x_b,y_a,y_b,nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE STIFFNESS AND LUMPED MASS MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA,MM,dmm] = get_mat_Q1_2d_xyl(node,elem,h,d,d_Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FINITE ELEMENT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fem = struct('node',node,'elem',elem,'AA',AA,'MM',MM,'dmm',dmm,...
             'nx',nx,'ny',ny,'h',h,'nno',nno,'dofs',dofs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIAL CONDITIONS OF THE STATE VARIABLES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1_0 = 100*ones(nno,1);
i1_0 = zeros(nno,1);
s2_0 = 50*ones(nno,1);
i2_0 = zeros(nno,1);
Z1_0 = zeros(nno,1);
for i=1:nno
   x=node(i,1);
   y=node(i,2);

%   s1_0(i)=1e+3*exp(-20*(x-1.5)^2-20*(y-0.2)^2);
   i1_0(i)=2e+1*exp(-100*(x-3.8)^2-100*(y-0.4)^2);
%   s2_0(i)=1e+3*exp(-20*(x-1.5)^2-20*(y-0.2)^2);
%   i2_0(i)=1e+1*exp(-100*(x-1.5)^2-100*(y-0.2)^2);
end
s1_0 = s1_0-i1_0;
uu0 = zeros(dofs,1);
uu0(1:5:end) = s1_0;
uu0(2:5:end) = i1_0;
uu0(3:5:end) = s2_0;
uu0(4:5:end) = i2_0;
uu0(5:5:end) = Z1_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE SYSTEM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muu = solve_sys(uu0,control,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_plot(muu,nx1,ny1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function M = M_fun(node)
%
% This functions returns the weed distribution  
%
% INPUT
%
% node      := node matrix, containing the x- and y-coordinates of 
%              the discretizations nodes 
%
% OUTPUT
%
% M         := weed distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = M_fun(node)

x=node(:,1);
y=node(:,2);

M=ones(size(x));

%for i=1:size(x,1)
%   if(x(i)>=1.5 && x(i)<=2.5)
%   if(x(i)>=2.5)
%     M(i)=0.4;
%   end
%end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ind_om = ind_om_fun(node)
%
% This functions returns the characteristic function of the weed distribution                             
%
% INPUT
%
% node      := node matrix, containing the x- and y-coordinates of
%              the discretizations nodes
%
% OUTPUT
%
% ind_om    := characteristic function of the control region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ind_om = ind_om_fun(node)

x=node(:,1);
y=node(:,2);

ind_om=zeros(size(x));

for i=1:size(x,1)
%   if(x(i)>=1.5 && x(i)<=2.5)
   if(x(i)>=2.5)
     ind_om(i)=1;
   end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function muu = solve_sys(uu0,control,param,fem,dt,tt)
%
% This function solves the system of PDEs 
%
% INPUT
%
% uu0      := initial condition for the state variables
% control  := control values
% param    := parameters
% fem      := finite element structure
% dt       := timestep size
% tt       := vector of timesteps
%
% OUTPUT
%
% muu      := spatio-temporal evolution of the state variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function muu = solve_sys(uu0,control,param,fem,dt,tt) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = fem.nx;
ny = fem.ny;
nx1 = nx+1;
ny1 = ny+1;
h = fem.h;
nno = nx1*ny1;
nel = nx*ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0;
k = 0; 
nt = length(tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIAL CONDITIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muu = zeros(fem.dofs,ceil(nt/100));
muu(:,1) = uu0;
uu = uu0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD ITERATION MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat = 1/dt*fem.MM+fem.AA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD M_FUN := function of weed distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = M_fun(fem.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD IND_OM_FUN := characteristic function of contro region 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind_om = ind_om_fun(fem.node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START TIME LOOP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = 1;
for k = 1:nt-1
  t = t+dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE REACTION TERM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  bbr = comp_reac(uu,control,param,fem,M,ind_om,t);  
  bb = 1/dt*fem.MM*uu+bbr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE LINEAR SYSTEM AND UPDATE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  uu = mat\bb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE STATE VARIABLES IN THE OUTPUT MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(mod(k,100)==1)
     kk = kk+1;
     muu(:,kk) = uu;
  end

end

clear bbr;
clear bb;  
clear mat;
clear M;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bbr = comp_reac(uu,control,param,fem,M,ind_om)
%
% This function computes the reaction term 
%
% INPUT
%
% uu       := vector of the state variables
% control  := control values
% param    := parameters
% fem      := finite element structure
% M        := weed distribution 
% ind_om   := characteristic function of the control region 
% t        := time
%
% OUTPUT
%
% bbr      := rection term 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bbr = comp_reac(uu,control,param,fem,M,ind_om,t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = param(2);
chi = param(3);
n = param(4);
q = param(5);
l = param(6);
C = param(7);
zeta = param(8);
b = param(9);
alpha = param(10);
mu = param(11);
p = param(12);
c = param(13);
r_2 = param(14);
d_Z = param(15);
delta_1 = param(16);
delta_2 = param(17);

K = 1;
beta = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT CONTROLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_21 = control(1);
gamma_1 = control(2);
gamma_23 = control(3);
gamma_22 = control(4);
gamma_3 = control(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT NUMBER OF NODES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nno = fem.nno;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT STATE VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s_1 = uu(1:5:end);
i_1 = uu(2:5:end);
s_2 = uu(3:5:end);
i_2 = uu(4:5:end);
Z_1 = uu(5:5:end);

C_I = s_1+i_1;
C_T = s_2+i_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE REACTION TERMS R_S1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rs_1 = zeros(nno,1);
k_i2 = rt1_fun(fem.node,fem.dmm,i_2);
rs_1 = r*C_I.*(M-chi*s_1)-n*s_1-p*delta_1*C_I.^2./(1+delta_2*C_I.^2).*...
       Z_1.*s_1./C_I-beta*s_1.*k_i2;
rs_1 = fem.dmm.*rs_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE REACTION TERMS R_I1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ri_1 = zeros(nno,1);
ri_1 = -n*i_1-p*delta_1*C_I.^2./(1+delta_2*C_I.^2).*Z_1.*i_1./C_I-...
       r*chi*i_1.*C_I+beta*s_1.*k_i2;
ri_1 = fem.dmm.*ri_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE REACTION TERMS R_S2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rs_2 = zeros(nno,1);
rs_2 = (q-l)*s_2-s_2.*C_T/C-(zeta*i_1+b*l*i_2).*s_2+alpha*i_2;
rs_2 = fem.dmm.*rs_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE REACTION TERMS R_I2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ri_2 = zeros(nno,1);
ri_2 = -(mu+l)*i_2-i_2.*C_T/C+(zeta*i_1+b*l*i_2).*s_2-alpha*i_2;
ri_2 = fem.dmm.*ri_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE REACTION TERMS R_Z1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rZ_1 = zeros(nno,1);
rZ_1 = r_2*Z_1.*(1-Z_1./(K+c*M))+delta_1*C_I.^2./(1+delta_2*C_I.^2).*Z_1;
if(mod(t,0.01)<0.003)
  rZ_1 = rZ_1+2e+5*ind_om;
  disp(['t = ',num2str(t),'  Z_1(10) = ',num2str(Z_1(10))]);
end
rZ_1 = fem.dmm.*rZ_1;

bbr = zeros(fem.dofs,1);
bbr(1:5:end) = rs_1;
bbr(2:5:end) = ri_1;
bbr(3:5:end) = rs_2;
bbr(4:5:end) = ri_2;
bbr(5:5:end) = rZ_1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bbr1 = rt1_fun(node,cc,i_2)
%
% This function computes the nonlocal term in the reaction terms of equations (6-7) 
%
% INPUT
%
% node   := node matrix, containing the x- and y-coordinates of 
%           the discretizations nodes
% cc     := diagonal mass matrix, obtained using mass lumping, in vector form
% i_2    := i_2 state variable
%
% OUTPUT
%
% bbr1   := nonlocal term 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bbr1=rt1_fun(node,cc,i_2)

nno=size(node,1);

bbr1=zeros(nno,1);

sigma=0.1;

%for i=1:nno
%   x1=node(i,1);
%   y1=node(i,2);
%   kk=zeros(nno,1);
%   for j=1:nno
%      x2=node(j,1);
%      y2=node(j,2);
%      kk(j)=exp((-(x1-x2)^2-(y1-y2)^2)/sigma);
%   end
%   bbr1(i)=dot(kk,cc.*i_2);
%end

bbr1 = i_2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [node,elem] = get_mesh_Q1_2d(x_a,x_b,y_a,y_b,nx,ny)
%
% This functions constructs the data structure of a 2D mesh for Q1 finite elements
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
% elem := elem (connectivity) matrix ((nx*ny)*4), containing for each element, 
%         the indices of the four vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
      node((j-1)*nx1+i,1) = x_a+(i-1)*hx;
      node((j-1)*nx1+i,2) = y_a+(j-1)*hy;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA,MM,cc] = get_mat_Q1_2d_xyl(node,elem,h,sigma,sigma_Z)
%
% This functions constructs the stiffness and mass matrices for 2D Q1 finite elements
%
% INPUT
%
% node   := node matrix, containing the x- and y-coordinates of 
%           the discretizations nodes
% elem   := elem (connectivity) matrix, containing for each element, 
%           the indices of the four vertices
% h      := mesh size (length of the generic sub-interval)
% sigma  := diffusion constant
%
% OUTPUT
%
% AA := block (4x4) stiffness matrix
% MM := block (4x4) diagonal mass matrix, obtained using mass lumping 
% cc := diagonal mass matrix, obtained using mass lumping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AA,MM,cc] = get_mat_Q1_2d_xyl(node,elem,h,sigma,sigma_Z) 

% nno := number of nodes
nno = size(node,1);

% nel := number of Q1 finite elements
nel = size(elem,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE OUTPUT MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA = sparse(5*nno,5*nno);
MM = sparse(5*nno,5*nno);
cc = zeros(nno,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE LOCAL MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA_loc,cc_loc] = get_mat_Q1_2d_loc(h);
MM_loc = diag(cc_loc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLING LOOP ON THE ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nel
   l = elem(i,:);

   kpm = l;
   cc(kpm) = cc(kpm)+cc_loc;

   kpm = 5*(l-1)+1;
   AA(kpm,kpm) = AA(kpm,kpm)+sigma*AA_loc;
   MM(kpm,kpm) = MM(kpm,kpm)+MM_loc;

   kpm = 5*(l-1)+2;
   AA(kpm,kpm) = AA(kpm,kpm)+sigma*AA_loc;
   MM(kpm,kpm) = MM(kpm,kpm)+MM_loc;

   kpm = 5*(l-1)+3;
   MM(kpm,kpm) = MM(kpm,kpm)+MM_loc;

   kpm = 5*(l-1)+4;
   MM(kpm,kpm) = MM(kpm,kpm)+MM_loc;

   kpm = 5*(l-1)+5;
   AA(kpm,kpm) = AA(kpm,kpm)+sigma_Z*AA_loc;
   MM(kpm,kpm) = MM(kpm,kpm)+MM_loc;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA_loc,cc_loc] = get_mat_Q1_2d_loc(h)
%
% This functions contsructs the local stiffness and mass matrices for 2D Q1 finite elements
%
% INPUT
%
% h      := mesh size (length of the generic sub-interval)
%
% OUTPUT
%
% AA_loc := local stiffness matrix
% cc_loc := local diagonal mass matrix, obtained using mass lumping, in vector form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AA_loc,cc_loc] = get_mat_Q1_2d_loc(h)

AA_loc = zeros(4,4);
cc_loc = zeros(4,1);

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

for i = 1:4
   cc_loc(i) = h^2/4;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_plot(muu,nx1,ny1)
%
% This functions makes the plot of the solution at three time instants 
%
% INPUT
%
% muu    := spatio-temporal evolution of the state variables
% nx1    := number of nodes on the x-direction
% ny1    := number of nodes on the y-direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_plot(muu,nx1,ny1)

figure(1)
i=101;
time=(i-1)*0.01;
s1=muu(1:5:end,i);
i1=muu(2:5:end,i);
s2=muu(3:5:end,i);
i2=muu(4:5:end,i);
plot_uu_stat_surf(s1,i1,s2,i2,nx1,ny1,time);

figure(2)
i=501;
time=(i-1)*0.01;
s1=muu(1:5:end,i);
i1=muu(2:5:end,i);
s2=muu(3:5:end,i);
i2=muu(4:5:end,i);
plot_uu_stat_surf(s1,i1,s2,i2,nx1,ny1,time);

figure(3)
i=1001;
time=(i-1)*0.01;
s1=muu(1:5:end,i);
i1=muu(2:5:end,i);
s2=muu(3:5:end,i);
i2=muu(4:5:end,i);
plot_uu_stat_surf(s1,i1,s2,i2,nx1,ny1,time);

figure(4)
i=1501;
time=(i-1)*0.01;
s1=muu(1:5:end,i);
i1=muu(2:5:end,i);
s2=muu(3:5:end,i);
i2=muu(4:5:end,i);
plot_uu_stat_surf(s1,i1,s2,i2,nx1,ny1,time);

figure(5)
i=2001;
time=(i-1)*0.01;
s1=muu(1:5:end,i);
i1=muu(2:5:end,i);
s2=muu(3:5:end,i);
i2=muu(4:5:end,i);
plot_uu_stat_surf(s1,i1,s2,i2,nx1,ny1,time);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plot_uu_stat_surf(uu_1,uu_2,uu_3,uu_4,nx1,ny1,time)
%
% This functions makes the plot of the solution at a given time instant 
%
% INPUT
%
% uu_1   := s_1 state variable
% uu_2   := i_1 state variable
% uu_3   := s_2 state variable
% uu_4   := i_2 state variable
% nx1    := number of nodes on the x-direction
% ny1    := number of nodes on the y-direction
% time   := time instant to plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_uu_stat_surf(uu_1,uu_2,uu_3,uu_4,nx1,ny1,time)

xx=linspace(0,400,nx1);
yy=linspace(0,80,ny1);

for j=1:ny1
   for i=1:nx1
      iind=(j-1)*nx1+i;
      mat_uu_3(j,i)=uu_3(iind);
      mat_uu_4(j,i)=uu_4(iind);
      mat_xx(j,i)=xx(i);
      mat_yy(j,i)=yy(j);
   end
end

h1=subplot(2,1,1);
surf(mat_xx,mat_yy,mat_uu_3)
axis([0 400 0 80 0 50])
set(gca,'XTick',[0:100:400],'YTick',[0:20:80],'ZTick',[0:10:50],'Fontsize',24)
xlabel('Km','Fontsize',20) 
ylabel('Km','Fontsize',20)
title(['s_2 (t = ',num2str(time),')'],'Fontsize',24)
view(6,16)

h2=subplot(2,1,2);
surf(mat_xx,mat_yy,mat_uu_4)
axis([0 400 0 80 0 15])
set(gca,'XTick',[0:100:400],'YTick',[0:20:80],'ZTick',[0:5:15],'Fontsize',24)
xlabel('Km','Fontsize',20)
ylabel('Km','Fontsize',20)
title(['i_2 (t = ',num2str(time),')'],'Fontsize',24)
view(6,16)

end
