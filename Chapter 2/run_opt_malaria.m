%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM RUN_OPT_MALARIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program solves the optimal control problem related to the spread of the
% malaria epidemic, paragraph 2.4,
% on the square domain [x_a,x_b] x [y_a,y_b]
% and the time interval [0,Tend].
% With the current set of parameters, it reproduces in particular
% the results reported in the second part of Table 2.2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1e-4;
d_22 = 0;
a_22 = 0.5;
eta = 10;
alpha_4 = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a = 0;
x_b = 1;
y_a = 0;
y_b = 1;
nx = 80;         % number of mesh elements in x-direction
ny = 80;         % number of mesh elements in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tend = 1;        % final time instant
dt = 0.05;       % time step size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONTROL VARIABLES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gam0_1 = 25.397; 
gam0_2 = 182.68; 
gam0_3 = 0.4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param(1) = d;
param(2) = d_22;
param(3) = a_22;
param(4) = eta;
param(5) = alpha_4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET CONTROL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma_1 = 400;
Gamma_2 = 400;
Gamma_3 = 0.80;
lim_gam_1 = 0;
lim_gam_2 = 0;
lim_gam_3 = 0;

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
[AA_3,dmm] = get_mat_Q1_2d(node,elem,h);

AA_1 = d*AA_3;
AA_2 = d_22*AA_3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FINITE ELEMENT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fem = struct('node',node,'elem',elem,'AA_1',AA_1,'AA_2',AA_2,'AA_3',AA_3,'dmm',dmm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIAL CONDITIONS FOR STATE VARIABLES u_1(x,t) and u_2(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu0_1 = zeros(nno,1);
uu0_2 = zeros(nno,1);
for i = 1:nno
   x = node(i,1);
   y = node(i,2);

   uu0_1(i) = 1e+3*exp(-4*(x-0.5)^2-4*(y-0.5)^2);
   uu0_2(i) = 1e+3*exp(-5*(x-0.5)^2-5*(y-0.5)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE CONTROLS gamma_1(t), gamma_2(t), gamma_3(t), phi(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gam_1 = gam0_1*ones(nt,1);
gam_2 = gam0_2*ones(nt,1);
gam_3 = gam0_3*ones(nt,1);

gams_1 = gam0_1;
gams_2 = gam0_2;
gams_3 = gam0_3;

phi = zeros(nno,1);
for i = 1:nno
    x = node(i,1);
    y = node(i,2);
    phi(i) = exp(-4*(x-0.5)^2-4*(y-0.5)^2)-0.75;
end
phi0 = phi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE FUNCTIONS HPhi, C, alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HPhi = HPhi_fun(phi);
C = C_fun(node);
alpha = alpha_fun(node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thres = 1e-3;
epsit = 1e-2;
maxiter = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE FORWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[muu_1,muu_2] = solve_pforw(uu0_1,uu0_2,gam_1,gam_2,gam_3,phi,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE BACKWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mqq_1,mqq_2] = solve_pback(muu_1,muu_2,gam_1,gam_2,gam_3,phi,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE COST FUNCTIONAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = comp_costJ(gam_1,gam_2,gam_3,phi,muu_1,muu_2,C,HPhi,alpha,tt,dmm,AA_3,param);

format shortG
disp([0,J,dot(dmm,muu_1(:,end)) dot(dmm,muu_2(:,end)),gam0_1,gam0_2,gam0_3])

errJ = 1;
i = 0;

while(errJ>1e-3 && i<maxiter)

   i=i+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE DERIVATIVES WITH RESPECT TO CONTROLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   F_1 = comp_F1(gam_1,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt);

   F_2 = comp_F2(gam_2,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt);

   F_3 = comp_F3(gam_3,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt);

   Ft = comp_Ft(gam_1,gam_2,gam_3,phi,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,alpha,fem,tt,param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE NEW CONTROLS AND CHECK DESCENT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   itera = 1;
   epsit_i = epsit;
   it_i = 0;

   while(itera && it_i<20)

      it_i = it_i+1;

      gam_1n = 0;
      gam_1n = gams_1-epsit_i*F_1;
      if(gam_1n<lim_gam_1)
        gam_1n = gams_1+epsit_i*F_1;
      end
      if(gam_1n>Gamma_1)
        gam_1n = Gamma_1-thres;
      end

      gam_2n = 0;
      gam_2n = gams_2-epsit_i*F_2;
      if(gam_2n<lim_gam_2)
        gam_2n = gams_2+epsit_i*F_2; 
      end
      if(gam_2n>Gamma_2)
        gam_2n = Gamma_2-thres;
      end

      gam_3n = 0;
      gam_3n = gams_3-epsit_i*F_3;
      if(gam_3n<lim_gam_3)
        gam_3n = gams_3+epsit_i*F_3;
      end
      if(gam_3n>Gamma_3)
        gam_3n = Gamma_3-thres;
      end

      gam_1 = gam_1n*ones(nt,1);
      gam_2 = gam_2n*ones(nt,1);
      gam_3 = gam_3n*ones(nt,1);

      phi_n = phi-epsit_i*Ft./fem.dmm;

      HPhi=HPhi_fun(phi_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE FORWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [muu_1,muu_2] = solve_pforw(uu0_1,uu0_2,gam_1,gam_2,gam_3,phi_n,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE BACKWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [mqq_1,mqq_2] = solve_pback(muu_1,muu_2,gam_1,gam_2,gam_3,phi_n,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE COST FUNCTIONAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Jn = comp_costJ(gam_1,gam_2,gam_3,phi_n,muu_1,muu_2,C,HPhi,alpha,tt,dmm,AA_3,param);

      if(Jn<=J)
        itera = 0;
      else
        epsit_i = epsit_i/10;
      end

   end

   errJ = abs(Jn-J)/abs(J);

   J = Jn;

   Jv(i) = J;

   epsit = epsit_i*100;

   mosq(i) = dot(dmm,muu_1(:,end));
   hum(i) = dot(dmm,muu_2(:,end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE CONTROLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   gams_1 = gam_1n;
   gams_2 = gam_2n;
   gams_3 = gam_3n;
   gam_1 = gams_1*ones(nt,1);
   gam_2 = gams_2*ones(nt,1);
   gam_3 = gams_3*ones(nt,1);
   phi = phi_n;

   disp([i,J,dot(dmm,muu_1(:,end)),dot(dmm,muu_2(:,end)),gams_1,gams_2,gams_3])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_uu(muu_1(:,end),muu_2(:,end),nx1,ny1);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function HPhi = HPhi_fun(phi)
%
% This function computes HPhi(phi) 
%
% INPUT
%
% phi  := distribution of control phi(x)
%
% OUTPUT
%
% HPhi   := function HPhi(phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HPhi = HPhi_fun(phi)

epsi = 1e-2;

HPhi = 1/2*(1+2/pi*atan(phi/epsi));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function C = C_fun(node)
%
% This function computes C(x) 
%
% INPUT
%
% node   := node matrix, containing the x- and y-coordinates of
%           the discretizations nodes
%
% OUTPUT
%
% C      := function C(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C=C_fun(node)

x=node(:,1);
y=node(:,2);

C=1e+4*exp(-5*(x-0.45).^2/30-5*(y-0.55).^2/30);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function alpha = alpha_fun(node)
%
% This function computes alpha(x)
%
% INPUT
%
% node   := node matrix, containing the x- and y-coordinates of
%           the discretizations nodes
%
% OUTPUT
%
% alpha  := function alpha(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = alpha_fun(node)

x = node(:,1);
y = node(:,2);

%alpha = 1e+3*x+1e+2*y;
alpha = 1e+2*ones(size(x));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function etv = eta_fun(node,eta)
%
% This function computes eta(x)
%
% INPUT
%
% node   := node matrix, containing the x- and y-coordinates of
%           the discretizations nodes
% eta    := parameter eta
%
% OUTPUT
%
% etv    := function eta(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etv = eta_fun(node,eta)

x = node(:,1);
y =node(:,2);

%etv = eta*(0.3*x+y)+0.08;
etv = eta*ones(size(x));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function F_1 = comp_F1(gam_1,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt)
%
% This function computes the derivative of the cost functional 
% with respect to control gamma_1
%
% INPUT
%
% gam_1 := control gamma_1
% dmm   := lumped mass matrix
% muu_1 := spatio-temporal evolution of u_1(x,t)
% muu_2 := spatio-temporal evolution of u_2(x,t)
% mqq_1 := spatio-temporal evolution of q_1(x,t)
% mqq_2 := spatio-temporal evolution of q_2(x,t)
% C     := function C
% HPhi  := function HPhi
% tt    := vector of timesteps
%
% OUTPUT
%
% F_1   := derivative of the cost functional with respect to control gamma_1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_1 = comp_F1(gam_1,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt)

Fv_1 = zeros(size(gam_1));

nt = size(muu_1,2);

for i = 1:nt
   iF_1 = HPhi.*(muu_1(:,i).*mqq_1(:,i)+C.*dzeta1(gam_1(i)));
   Fv_1(i) = dot(dmm,iF_1);
end

F_1 = trapz(tt,Fv_1);

clear iF_1
clear Fv_1

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function F_2 = comp_F2(gam_2,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt)
%
% This function computes the derivative of the cost functional 
% with respect to control gamma_2
%
% INPUT
%
% gam_2 := control gamma_2
% dmm   := lumped mass matrix
% muu_1 := spatio-temporal evolution of u_1(x,t)
% muu_2 := spatio-temporal evolution of u_2(x,t)
% mqq_1 := spatio-temporal evolution of q_1(x,t)
% mqq_2 := spatio-temporal evolution of q_2(x,t)
% C     := function C
% HPhi  := function HPhi
% tt    := vector of timesteps
%
% OUTPUT
%
% F_2   := derivative of the cost functional with respect to control gamma_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_2 = comp_F2(gam_2,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt)

Fv_2 = zeros(size(gam_2));

nt = size(muu_1,2);

for i = 1:nt
   iF_2 = HPhi.*(muu_2(:,i).*mqq_2(:,i)+muu_2(:,i).*dzetat2(gam_2(i)*muu_2(:,i)));
   Fv_2(i) = dot(dmm,iF_2);
end

F_2 = trapz(tt,Fv_2);

clear iF_2
clear Fv_2

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function F_3 = comp_F3(gam_3,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt)
%
% This function computes the derivative of the cost functional
% with respect to control gamma_3
%
% INPUT
%
% gam_3 := control gamma_3
% dmm   := lumped mass matrix
% muu_1 := spatio-temporal evolution of u_1(x,t)
% muu_2 := spatio-temporal evolution of u_2(x,t)
% mqq_1 := spatio-temporal evolution of q_1(x,t)
% mqq_2 := spatio-temporal evolution of q_2(x,t)
% C     := function C
% HPhi  := function HPhi
% tt    := vector of timesteps
%
% OUTPUT
%
% F_3   := derivative of the cost functional with respect to control gamma_3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F_3 = comp_F3(gam_3,dmm,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,tt)

F_3 = zeros(size(gam_3));

nt = size(muu_1,2);

for i = 1:nt
   iF_3 = HPhi.*((C-muu_2(:,i)).*muu_1(:,i)./C.*mqq_2(:,i)+C.*dzeta3(gam_3(i)));
   Fv_3(i) = dot(dmm,iF_3);
end

F_3 = trapz(tt,Fv_3);

clear iF_3
clear Fv_3

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function J = comp_costJ(gam_1,gam_2,gam_3,phi,muu_1,muu_2,C,HPhi,alpha,tt,dmm,AA_3,param)
%
% This function computes the cost functional
%
% INPUT
%
% gam_1 := control gamma_1
% gam_2 := control gamma_2
% gam_3 := control gamma_3
% phi   := control phi 
% muu_1 := spatio-temporal evolution of u_1(x,t)
% muu_2 := spatio-temporal evolution of u_2(x,t)
% C     := function C
% HPhi  := function HPhi
% alpha := function alpha 
% tt    := vector of timesteps
% dmm   := lumped mass matrix
% AA_3  := stiffness matrix 
% param := parameters
%
% OUTPUT
%
% J     := cost functional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = comp_costJ(gam_1,gam_2,gam_3,phi,muu_1,muu_2,C,HPhi,alpha,tt,dmm,AA_3,param)

alpha_4 = param(5);;

nt = size(muu_1,2);

for i = 1:nt
   iJJ = C.*zeta1(gam_1(i)).*HPhi+...
         zeta2(gam_2(i)*muu_2(:,i)).*gam_2(i).*muu_2(:,i).*HPhi+...
         C.*zeta3(gam_3(i)).*HPhi+...
         alpha_4*muu_2(:,i);
   iJ(i) = dot(dmm,iJJ);
end
J = trapz(tt,iJ);

iJJ = alpha.*HPhi;
J = J+dot(dmm,iJJ);

epsilon = 1e-2;
beta = 1e-2;

J = J+beta*epsilon*HPhi'*AA_3*HPhi;

iJJ = HPhi.^2.*(1-HPhi).^2;
J = J+beta/epsilon*dot(dmm,iJJ);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [muu_1,muu_2] = solve_pforw(uu0_1,uu0_2,gam_1,gam_2,gam_3,phi,param,fem,dt,tt)
%
% This function solves the forward problem 
%
% INPUT
%
% uu0_1 := initial distribution of u_1(x,t)
% uu0_2 := initial distribution of u_2(x,t)
% gam_1 := control gamma_1
% gam_2 := control gamma_2
% gam_3 := control gamma_3
% phi   := control phi
% param := parameters
% fem   := finite element structure
% dt    := timestep size
% tt    := vector of timesteps
%
% OUTPUT
%
% muu_1 := spatio-temporal evolution of u_1(x,t)
% muu_2 := spatio-temporal evolution of u_2(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [muu_1,muu_2] = solve_pforw(uu0_1,uu0_2,gam_1,gam_2,gam_3,phi,param,fem,dt,tt)

% set parameters
d = param(1);
d_22 = param(2);
a_22 = param(3);
eta = param(4);

% set mesh parameters
nno = size(fem.node,1);

% set time loop parameters
t = 0;
k = 0; 
dt = tt(2)-tt(1);
nt = length(tt);

% set initial conditions
muu_1 = zeros(nno,nt);
muu_2 = zeros(nno,nt);
muu_1(:,1) = uu0_1;
muu_2(:,1) = uu0_2;
uu_1 = uu0_1;
uu_2 = uu0_2;

% build HPhi fun
HPhi = HPhi_fun(phi);

Mass = spdiags(fem.dmm,0,nno,nno);
HHPhi = spdiags(fem.dmm.*HPhi,0,nno,nno);

% build C fun
C = C_fun(fem.node);

% build eta fun
etv = eta_fun(fem.node,eta);

Metv = spdiags(fem.dmm.*etv,0,nno,nno);

% time loop start
for k = 1:nt-1

   t = t+dt;

   % build iteration matrix
   mat = [1/dt*Mass+Metv+fem.AA_1+gam_1(k+1)*HHPhi sparse(nno,nno);
          sparse(nno,nno)                          (1/dt+a_22)*Mass+fem.AA_2+gam_2(k+1)*HHPhi];


   rt_1 = rt1_fun(fem.node,fem.dmm,uu_2);
   bb_1 = fem.dmm.*(1/dt*uu_1+rt_1);

   G = G_fun2(C,uu_1,uu_2);
   bb_2 = fem.dmm.*(1/dt*uu_2+(1-gam_3(k+1).*HPhi).*G);
 
   bb = [bb_1;bb_2];

   uu = mat\bb;

% update
   uu_1 = uu(1:nno);
   uu_2 = uu(nno+1:end);

   muu_1(:,k+1) = uu_1;
   muu_2(:,k+1) = uu_2;

end

clear uu_1;
clear uu_2;  
clear rt_1;
clear HPhi;
clear HHPhi;
clear Mass;
clear mat;
clear bb_1;
clear bb_2;
clear bb;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bbr = rt1_fun(node,dmm,uu_2)
%
% This function computes the reaction term in the first equaion of the 
% forward problem 
%
% INPUT
%
% node   := node matrix, containing the x- and y-coordinates of
%           the discretizations nodes
% dmm    := lumped mass matrix
% uu_2   := spatial distribution of u_2(x,t) at current time step
%
% OUTPUT
%
% bbr    := reaction term evaluated in each mesh node 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bbr = rt1_fun(node,dmm,uu_2)

nno = size(node,1);

bbr = zeros(nno,1);

sigma = 0.1;

for i =1:nno
   x1 = node(i,1);
   y1 = node(i,2);
   kk = zeros(nno,1);
   for j = 1:nno
      x2 = node(j,1);
      y2 = node(j,2);
      kk(j) = exp((-(x1-x2)^2-(y1-y2)^2)/sigma);
   end
   bbr(i) = dot(kk,dmm.*uu_2);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [mqq_1,mqq_2] = solve_pback(muu_1,muu_2,gam_1,gam_2,gam_3,phi,param,fem,dt,tt)
%
% This function solves the backward problem 
%
% INPUT
%
% muu_1 := spatio-temporal evolution of u_1(x,t)
% muu_2 := spatio-temporal evolution of u_2(x,t)
% gam_1 := control gamma_1
% gam_2 := control gamma_2
% gam_3 := control gamma_3
% phi   := control phi
% param := parameters
% fem   := finite element structure
% dt    := timestep size
% tt    := vector of timesteps
%
% OUTPUT
%
% mqq_1 := spatio-temporal evolution of q_1(x,t)
% mqq_2 := spatio-temporal evolution of q_2(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mqq_1,mqq_2] = solve_pback(muu_1,muu_2,gam_1,gam_2,gam_3,phi,param,fem,dt,tt)

% set parameters
d = param(1);
d_22 = param(2);
a_22 = param(3);
eta = param(4);
alpha_4 = param(5);

% set mesh parameters
nno = size(fem.node,1);

% set time loop parameters
t = 0;
k = 0;
dt = tt(2)-tt(1);
nt = length(tt);

% set initial conditions
mqq_1 = zeros(nno,nt);
mqq_2 = zeros(nno,nt);
qq_1 = mqq_1(:,end);
qq_2 = muu_2(:,end);

% build HPhi fun
HPhi = HPhi_fun(phi);
dHPhi = dHPhi_fun(phi);

% build C fun
C = C_fun(fem.node);
 
Mass = spdiags(fem.dmm,0,nno,nno);
HHPhi = spdiags(fem.dmm.*HPhi,0,nno,nno);
dHHPhi = spdiags(fem.dmm.*dHPhi,0,nno,nno);

% build eta fun
etv = eta_fun(fem.node,eta);

Metv = spdiags(fem.dmm.*etv,0,nno,nno);

% time loop start
for k = 1:nt-1

   t = tt(nt-k);
 
   % build iteration matrix
   mat = [1/dt*Mass+Metv+fem.AA_1+gam_1(nt-k)*dHHPhi   sparse(nno,nno);
          sparse(nno,nno)                              (1/dt+a_22)*Mass+fem.AA_2+gam_2(nt-k)*HHPhi];

   uu_1 = muu_1(:,nt-k);
   uu_2 = muu_2(:,nt-k); 

   bb_1 = fem.dmm.*(1/dt*qq_1+(1-gam_3(nt-k)*HPhi).*(C-uu_2)./C.*qq_2).*dhz_fun(uu_1./C); 

   rt_2 = rt1_fun(fem.node,fem.dmm,qq_1);
   rt_2 = -rt_2+(1-gam_3(nt-k)*HPhi).*hz_fun(uu_1./C).*qq_2;
   rt_2 = rt_2+dzetat2(gam_2(nt-k)*uu_2).*gam_2(nt-k).*HPhi+alpha_4;

   bb_2 = fem.dmm.*(1/dt*qq_2-rt_2); 
 
   bb = [bb_1;bb_2];

   qq = mat\bb;

% update
   qq_1 = qq(1:nno);
   qq_2 = qq(nno+1:end);

   mqq_1(:,nt-k) = qq_1;
   mqq_2(:,nt-k) = qq_2;

end

clear uu_1;
clear uu_2;  
clear qq_1;
clear qq_2;
clear rt_2;
clear HPhi;
clear HHPhi;
clear Mass;
clear mat;
clear bb_1;
clear bb_2;
clear bb;
clear dHPhi
clear dHHPhi;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Ft = comp_Ft(gam_1,gam_2,gam_3,phi,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,alpha,fem,tt,param)
%
% This function computes the derivative of the cost functional
% with respect to control phi
%
% INPUT
%
% gam_1 := control gamma_1
% gam_2 := control gamma_2
% gam_3 := control gamma_3
% phi   := control phi
% muu_1 := spatio-temporal evolution of u_1(x,t)
% muu_2 := spatio-temporal evolution of u_2(x,t)
% mqq_1 := spatio-temporal evolution of q_1(x,t)
% mqq_2 := spatio-temporal evolution of q_2(x,t)
% C     := function C
% HPhi  := function HPhi
% alpha := function alpha
% fem   := fem structure
% tt    := vector of timesteps
% param := parameters
%
% OUTPUT
%
% Ft    := derivative of the cost functional with respect to control phi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ft = comp_Ft(gam_1,gam_2,gam_3,phi,muu_1,muu_2,mqq_1,mqq_2,C,HPhi,alpha,fem,tt,param)

alpha_4 = param(5);

nno = size(muu_1,1);
nt = size(muu_1,2);

% Compute Ft_1
Ft_1 = zeros(nno,1);

for j = 1:nno
   iF = zeros(nt,1);
   for i = 1:nt
      iF(i) = C(j)*zeta1(gam_1(i)) + zetat2(gam_2(i)*muu_2(j,i)) + C(j)*zeta3(gam_3(i))+...
              gam_1(i)*muu_1(j,i)*mqq_1(j,i) + gam_2(i)*muu_2(j,i)*mqq_2(j,i)+...
	      gam_3(i)*(C(j)-muu_2(j,i))*hz_fun(muu_1(j,i)/C(j))*mqq_2(j,i) + alpha(j);
   end
   Ft_1(j) = trapz(tt,iF);
end

% Compute Ft_2
epsilon = 1e-2;
beta = 1e-2;
sigma = 2*beta*epsilon;

Aphi = comp_Aphi(fem.elem,nno,sigma,phi);

Ft_2 = Aphi*phi;

% Compute Ft_3
Ft_3 = -beta/epsilon*(4*HPhi.^3-6*HPhi.^2+2*HPhi);

Ft = Ft_2+fem.dmm.*(Ft_1+Ft_3);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Aphi = comp_Aphi(elem,nno,sigma,phi)
%
% This function computes the stiffness matrix associated with the 
% nonlinear diffusion term in the derivative of the cost functional
% with respect to control phi
%
% INPUT
%
% elem   := elem (connectivity) matrix, containing for each
%           element, the indices of the four vertices
% nno    := number of nodes
% sigma  := diffusion coefficient
% phi    := control phi
%
% OUTPUT
%
% Aphi   := stiffness matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aphi = comp_Aphi(elem,nno,sigma,phi) 

nel = size(elem,1);

Aphi = sparse(nno,nno);

[AA_loc,cc_loc] = get_mat_Q1_2d_loc(1);
AA_loc = sigma*AA_loc;

for i = 1:nel
   l1 = elem(i,1);
   l2 = elem(i,2);
   l3 = elem(i,3);
   l4 = elem(i,4);

   kpm(1) = l1;
   kpm(2) = l2;
   kpm(3) = l3;
   kpm(4) = l4;	 

   phim=mean(phi(kpm));
   dHPhi=dHPhi_fun(phim);

   Aphi(kpm,kpm) = Aphi(kpm,kpm)+dHPhi*AA_loc;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUXILIARY FUNCTIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function G = G_fun2(C,uu_1,uu_2)
function G = G_fun2(C,uu_1,uu_2)

hz = hz_fun(uu_1./C);
G = (C-uu_2).*hz;

end


% function hz = hz_fun(z)
function hz = hz_fun(z)

hz = 4*z;

end


% function dhz = dhz_fun(z)
function dhz = dhz_fun(z)

dhz = 4;

end


% function dHPhi = dHPhi_fun(phi)
function dHPhi = dHPhi_fun(phi)

epsi = 1e-2;

dHPhi = 1/(epsi*pi)*1./(1+(phi/epsi).^2);

end


% function z1 = zeta1(s)
function z1 = zeta1(s)

Gamma_1 = 500;
b_1 = 1;
alpha_1 = 1e-1;

z1 = s./((Gamma_1-s).*b_1);
z1 = alpha_1*z1;

end


% function dz1 = dzeta1(s)
function dz1 = dzeta1(s)

Gamma_1 = 500;
b_1 = 1;
alpha_1 = 1e-1;

dz1 = Gamma_1./((Gamma_1-s).^2.*b_1);
dz1 = alpha_1*dz1;

end


% function z2 = zeta2(s)
function z2 = zeta2(s)

Gamma_2 = 500;
a_2 = 0.02;
c_2 = 600;
alpha_2 = 1e-2;

z2 = alpha_2*(Gamma_2*a_2*s+c_2)./(1+a_2*s);

end


% function dz2 = dzeta2(s)
function dz2 = dzeta2(s)

Gamma_2 = 500;
a_2 = 0.02;
c_2 = 600;
alpha_2 = 1e-2;

dz2 = alpha_2*(Gamma_2-c_2)*a_2./(1+a_2*s).^2;

end


% function zt2 = zetat2(s)
function zt2 = zetat2(s)

zt2 = zeta2(s).*s;

end


% function dzt2 = dzetat2(s)
function dzt2 = dzetat2(s)

dzt2 = dzeta2(s).*s+zeta2(s);

end


% function z3 = zeta3(s)
function z3 = zeta3(s)

Gamma_3 = 0.85;
b_3 = 1;
alpha_3 = 1e-2;

z3 = s./((Gamma_3-s).*b_3);
z3 = alpha_3*z3;

end


% function dz3 = dzeta3(s)
function dz3 = dzeta3(s)

Gamma_3 = 0.85;
b_3 = 1;
alpha_3 = 1e-2;

dz3 = Gamma_3./((Gamma_3-s).^2.*b_3);
dz3 = alpha_3*dz3;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function plot_uu(uu_1,uu_2,nx1,ny1) 
%
% This function plots the distributions of u_1(x,t) and u_2(x,t)
% at the current time step 
%
% INPUT
%
% uu_1 := spatial distribution of u_1(x,t)
% uu_2 := spatial distribution of u_2(x,t)
% nx1  := number of nodes along the x-direction
% ny1  := number of nodes along the y-direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_uu(uu_1,uu_2,nx1,ny1)

x=linspace(0,1,nx1);
y=linspace(0,1,ny1);

for j=1:ny1
   for i=1:nx1
      iind=(j-1)*nx1+i;
      mat_uu_1(j,i)=uu_1(iind);
      mat_uu_2(j,i)=uu_2(iind);
      mat_x(j,i)=x(i);
      mat_y(j,i)=y(j);
   end
end

h1=subplot(1,2,1);
surf(mat_x,mat_y,mat_uu_1)
axis([0 1 0 1 0 140])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1],'ZTick',[0:40:140],'Fontsize',30)
set(h1,'CLim',[0 140])
title('infected mosquitos','Fontsize',30)
colorbar('Fontsize',30)

h1=subplot(1,2,2);
surf(mat_x,mat_y,mat_uu_2)
axis([0 1 0 1 0 1400])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1],'ZTick',[0:400:1400],'Fontsize',30)
set(h1,'CLim',[0 1400])
title('infected humans','Fontsize',30)
colorbar('Fontsize',30)

end
