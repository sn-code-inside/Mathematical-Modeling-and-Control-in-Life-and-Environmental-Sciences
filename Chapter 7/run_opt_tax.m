%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM RUN_OPT_TAX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program solves the optimal control problem for the economic growth model 
% with pollution diffusion and environmental taxation,
% according to the conceptual algorithm reported in section 7.4.1,
% on the one-dimensional domain [x_a,x_b] and the time interval [0,Tend].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_1 = 1e-3;
d_2 = 1e-3;
delta_1 = 0.05;
delta_2 = 0.03;
A = 1;
s = 0.6;
theta = 2;
alpha_1 = 0.7;
alpha_2 = 1;
gamma = 4;
phi = 0.3;
lambda_0 = 1;
lambda_1 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a = -1;
x_b = 1;
nx = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tend = 10;
dt = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MODEL PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param(1) = d_1;
param(2) = d_2;
param(3) = delta_1;
param(4) = delta_2;
param(5) = A;
param(6) = s;
param(7) = theta;
param(8) = alpha_1;
param(9) = alpha_2;
param(10) = gamma;
param(11) = phi;
param(12) = lambda_0;
param(13) = lambda_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MESH PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1 = nx+1;
h = (x_b-x_a)/nx;
nno = nx1;
nel = nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET TIME PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = [0:dt:Tend];
nt = length(tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NODE AND ELEM MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[node,elem] = get_mesh_P1_1d(x_a,x_b,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE STIFFNESS AND LUMPED MASS MATRICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA,dmm] = get_mat_P1_1d(node,elem,h);

AA_1 = d_1*AA;
AA_2 = d_2*AA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FINITE ELEMENT STRUCTURE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fem = struct('node',node,'elem',elem,'AA_1',AA_1,'AA_2',AA_2,'dmm',dmm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIAL CONDITIONS FOR STATE VARIABLES k(x,t) and p(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk0 = zeros(nno,1);
pp0 = zeros(nno,1);
for i=1:nno
   x = node(i,1);

   kk0(i)=exp(-(x-0.5)^2);
   pp0(i)=3*exp(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE CONTROLS c(x,t) and tau(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcc = zeros(nno,nt);
mtau = zeros(nno,nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OPTIMIZATION ALGORITHM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err = 1;
op_it = 0;
epsi_0 = 0.001;

while(err>1e-3 && op_it<10000) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE FORWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('solve the forward problem')
  [mkk,mpp] = solve_pforw(kk0,pp0,mcc,mtau,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE BACKWARD PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('solve the backward problem')
  [mllk,mllp] = solve_pback(mkk,mpp,mcc,mtau,param,fem,dt,tt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE BETA_STAR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mbeta_st = zeros(nno,nt);

  for i = 1:nt
     cc = mcc(:,i);
     tau = mtau(:,i);
     llk = mllk(:,i);
     llp = mllp(:,i);
     kk = mkk(:,i);
     pp = mpp(:,i);

     mbeta_st(:,i) = comp_beta_st(llk,llp,kk,pp,cc,tau,param,dmm);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE COST FUNCTIONAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  I = comp_I(mkk,mpp,mcc,mtau,param,tt,dmm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE CONTROLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mv_1 = zeros(nno,nt);
  mv_2 = zeros(nno,nt);

  for i = 1:nno
     for j = 1:nt
        mv_1(i,j) = sign(1-mllk(i,j));
        mv_2(i,j) = -sign(mbeta_st(i,j));
     end
  end

  mcc_ = mcc+epsi_0*mv_1;
  mtau_ = mtau+epsi_0*mv_2;

  xi = linspace(0,1,100);
  for i = 1:length(xi)
     mcc_xi = xi(i)*mcc+(1-xi(i))*mcc_;
     mtau_xi = xi(i)*mtau+(1-xi(i))*mtau_;
     Ixi(i) = comp_I(mkk,mpp,mcc_xi,mtau_xi,param,tt,dmm);
  end

  [max_Ixi,imax_Ixi] = max(Ixi);
  xi_0 = xi(imax_Ixi);

  mcc_ = xi_0*mcc+(1-xi_0)*mcc_;
  mtau_ = xi_0*mtau+(1-xi_0)*mtau_;  

  for j = 1:nno
     for i = 1:nt
        mtau_(j,i) = max(mtau_(j,i),0);
        mtau_(j,i) = min(mtau_(j,i),0.4);
        if(1-mllk(j,i)<0)
          mcc_(j,i) = 0;
        else
          mcc_(j,i) = 1-s-mtau_(j,i);
        end
     end
  end

  err = max(norm(mcc_-mcc),norm(mtau_-mtau))/max(norm(mcc),norm(mtau))
  op_it = op_it+1

  mcc = mcc_;
  mtau = mtau_;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:101
   for j=1:201
      mat_x(i,j)=node(i);
      mat_t(i,j)=j*dt;
   end
end

figure(1)
surf(mat_x',mat_t',mkk')
axis([-1 1 0 10 0 1])
set(gca,'XTick',[-1:0.5:1],'YTick',[0:5:10],'ZTick',[0:0.5:1],'Fontsize',30)
xlabel('x','Fontsize',30)
ylabel('t','Fontsize',30)
title('k(x,t)','Fontsize',30)

figure(2)
surf(mat_x',mat_t',mpp')
axis([-1 1 0 10 0 9])
set(gca,'XTick',[-1:0.5:1],'YTick',[0:5:10],'ZTick',[0:3:9],'Fontsize',30)
xlabel('x','Fontsize',30)
ylabel('t','Fontsize',30)
title('p(x,t)','Fontsize',30)

figure(3)
surf(mat_x',mat_t',mcc')
axis([-1 1 0 10 0.397 0.4])
set(gca,'XTick',[-1:0.5:1],'YTick',[0:5:10],'ZTick',[0.397:0.001:0.4],'Fontsize',30)
xlabel('x','Fontsize',30)
ylabel('t','Fontsize',30)
title('c(x,t)','Fontsize',30)

figure(4)
surf(mat_x',mat_t',mtau')
axis([-1 1 0 10 0 0.003])
set(gca,'XTick',[-1:0.5:1],'YTick',[0:5:10],'ZTick',[0:0.001:0.003],'Fontsize',30)
xlabel('x','Fontsize',30)
ylabel('t','Fontsize',30)
title('\tau(x,t)','Fontsize',30)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [mkk,mpp] = solve_pforw(kk0,pp0,mcc,mtau,param,fem,dt,tt)
%
% This functions solves the forward problem, 
% i.e. the system of two equations for k(x,t) and p(x,t)
%
% INPUT
%
% kk0   := initial condition for the state variable k(x,t) 
% pp0   := initial condition for the state variable p(x,t) 
% mcc   := spatio-temporal evolution of control c(x,t) 
% mtau  := spatio-temporal evolution of control tau(x,t) 
% param := parameters
% fem   := finite element structure
% dt    := timestep size
% tt    := vector of timesteps
%
% OUTPUT
%
% mkk   := spatio-temporal evolution of the state variable k(x,t)
% mpp   := spatio-temporal evolution of the state variable p(x,t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mkk,mpp] = solve_pforw(kk0,pp0,mcc,mtau,param,fem,dt,tt)

% get the number of nodes 
nno = size(kk0,1);

% set time loop parameters
nt = length(tt);

% set initial conditions
mkk = zeros(nno,nt);
mpp = zeros(nno,nt);
mkk(:,1) = kk0;
mpp(:,1) = pp0;

% assemble iteration matrix
Mass = spdiags(fem.dmm,0,nno,nno);
mat = [1/dt*Mass+fem.AA_1   sparse(nno,nno);
       sparse(nno,nno)      1/dt*Mass+fem.AA_2];

% time loop start
for k=1:nt-1
   t = tt(k+1);

   % compute and integrate reaction term of equation for k(x,t)
   rt_1 = rt1_fun(mkk(:,k),mpp(:,k),mcc(:,k),mtau(:,k),param);
   bb_1 = fem.dmm.*(1/dt*mkk(:,k)+rt_1);

   % compute and integrate reaction term of equation for p(x,t)
   rt_2 = rt2_fun(mkk(:,k),mpp(:,k),mcc(:,k),mtau(:,k),param,fem.dmm);
   bb_2 = fem.dmm.*(1/dt*mpp(:,k)+rt_2);
 
   bb = [bb_1;bb_2];

   % solve the linear system
   uu=mat\bb;

   % update
   mkk(:,k+1) = uu(1:nno);
   mpp(:,k+1) = uu(nno+1:end);
end

clear uu;
clear rt_1;
clear rt_2;
clear Mass;
clear mat;
clear bb_1;
clear bb_2;
clear bb;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rt_1=rt1_fun(kk,pp,cc,tau,param); 
%
% This functions constructs the reaction term of the equation for k(x,t) 
%
% INPUT
%
% kk    := current state variable k(x,t) distribution
% pp    := current state variable p(x,t) distribution
% cc    := current control c(x,t) distribution
% tau   := current control tau(x,t) distribution 
% param := parameters
%
% OUTPUT
%
% rt_1  := reaction term of the first equation evaluated in each mesh node 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt_1=rt1_fun(kk,pp,cc,tau,param); 

% get the number of nodes
nno = size(kk,1);

% initialize the output
rt_1 = zeros(nno,1);

% extract useful parameters
delta_1 = param(3);
A = param(5);
alpha_1 = param(8);
alpha_2 = param(9);
gamma = param(10);

% compute the reaction term
ff = f_fun((1-tau).*kk,alpha_1,alpha_2,gamma);
rt_1 = A./(1+pp.^2).*ff-delta_1*kk-cc.*kk;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rt_2 = rt2_fun(kk,pp,cc,tau,param,dmm)
%
% This functions constructs the reaction term of the equation for p(x,t)
%
% INPUT
%
% kk    := current state variable k(x,t) distribution
% pp    := current state variable p(x,t) distribution
% cc    := current control c(x,t) distribution
% tau   := current control tau(x,t) distribution
% param := parameters
% dmm   := diagonal mass matrix, obtained using mass lumping 
%
% OUTPUT
%
% rt_1  := reaction term of the first equation evaluated in each mesh node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt_2 = rt2_fun(kk,pp,cc,tau,param,dmm)

% get the number of nodes
nno = size(kk,1);

% initialize the output
rt_2 = zeros(nno,1);

% extract useful parameters
delta_2 = param(4);
theta = param(7);
alpha_1 = param(8);
alpha_2 = param(9);
gamma = param(10);
phi = param(11);

% compute the non-local term
ff = f_fun((1-tau).*kk,alpha_1,alpha_2,gamma);
bb_nl = dot(ff,dmm)*phi*ones(nno,1);

% compute the reaction term
rt_2 = theta*bb_nl-delta_2*pp;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ff = f_fun(rr,alpha_1,alpha_2,gamma)
%
% This function evaluates the production function f
%
% INPUT
%
% rr      := amount of physical capital 
% alpha_1 := positive parameter 
% alpha_2 := positive parameter 
% gamma   := positive parameter 
%
% OUTPUT
%
% ff      := gross domestic production
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ff = f_fun(rr,alpha_1,alpha_2,gamma) 

ff = alpha_1*rr.^gamma./(1+alpha_2*rr.^gamma);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function dff = df_fun(rr,alpha_1,alpha_2,gamma)
%
% This function evaluates the derivative of the production function f
%
% INPUT
%
% rr      := amount of physical capital
% alpha_1 := positive parameter
% alpha_2 := positive parameter
% gamma   := positive parameter
%
% OUTPUT
%
% dff      := derivative of the gross domestic production
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dff = df_fun(rr,alpha_1,alpha_2,gamma)

dff = (alpha_1*gamma*rr.^(gamma-1).*(1+alpha_2*rr.^gamma)-alpha_1*rr.^gamma.*...
      alpha_2.*gamma.*rr.^(gamma-1))./(1+alpha_2*rr.^gamma).^2;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [mllk,mllp] = solve_pback(mkk,mpp,mcc,mtau,param,fem,dt,tt)
%
% This functions solves the backward problem, i.e. the system of two equations 
% for lambda_k(x,t) and lambda_p(x,t)
%
% INPUT
%
% mkk   := spatio-temporal evolution of the state variable k(x,t)
% mpp   := spatio-temporal evolution of the state variable p(x,t)
% mcc   := spatio-temporal evolution of control c(x,t)
% mtau  := spatio-temporal evolution of control tau(x,t)
% param := parameters
% fem   := finite element structure
% dt    := timestep size
% tt    := vector of timesteps
%
% OUTPUT
%
% mllk  := spatio-temporal evolution of lambda_k(x,t)
% mllp  := spatio-temporal evolution of lambda_p(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mllk,mllp] = solve_pback(mkk,mpp,mcc,mtau,param,fem,dt,tt)

% get the number of nodes
nno = size(mkk,1);

% set time loop parameters
nt = length(tt);

% set initial conditions
mllk = zeros(nno,nt);
mllp = zeros(nno,nt);

% assemble iteration matrix
Mass = spdiags(fem.dmm,0,nno,nno);
mat = [1/dt*Mass+fem.AA_1   sparse(nno,nno);
       sparse(nno,nno)      1/dt*Mass+fem.AA_2];

% time loop start
for k=1:nt-1
   t = tt(end-k);

   % compute and integrate reaction term of equation for lambda_k(x,t)
   rt_1 = rtb1_fun(mllk(:,nt-k+1),mllp(:,nt-k+1),mkk(:,nt-k),mpp(:,nt-k),...
          mcc(:,nt-k),mtau(:,nt-k),param,fem.dmm);
   bb_1 = fem.dmm.*(1/dt*mllk(:,nt-k+1)-rt_1);

   % compute and integrate reaction term of equation for lambda_p(x,t)
   rt_2 = rtb2_fun(mllk(:,nt-k+1),mllp(:,nt-k+1),mkk(:,nt-k),mpp(:,nt-k),...
          mcc(:,nt-k),mtau(:,nt-k),param);
   bb_2 = fem.dmm.*(1/dt*mllp(:,nt-k+1)-rt_2);

   bb = [bb_1;bb_2];

   % solve the linear system
   uu=mat\bb;

   % update
   mllk(:,nt-k) = uu(1:nno);
   mllp(:,nt-k) = uu(nno+1:end);
end

clear uu;
clear rt_1;
clear rt_2;
clear Mass;
clear mat;
clear bb_1;
clear bb_2;
clear bb;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rt_1 = rtb1_fun(llk,llp,kk,pp,cc,tau,param,dmm)
%
% This functions constructs the reaction term of the equation for lambda_k(x,t)
%
% INPUT
%
% llk   := current lambda_k(x,t) distribution
% llp   := current lambda_p(x,t) distribution
% kk    := current state variable k(x,t) distribution
% pp    := current state variable p(x,t) distribution
% cc    := current control c(x,t) distribution
% tau   := current control tau(x,t) distribution
% param := parameters
% dmm   := diagonal mass matrix, obtained using mass lumping
%
% OUTPUT
%
% rt_1  := reaction term of the first equation evaluated in each mesh node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt_1 = rtb1_fun(llk,llp,kk,pp,cc,tau,param,dmm)

% get the number of nodes
nno = size(kk,1);

% initialize the output
rt_1 = zeros(nno,1);

% extract useful parameters
delta_1 = param(3);
A = param(5);
theta = param(7);
alpha_1 = param(8);
alpha_2 = param(9);
gamma = param(10);
phi = param(11);
lambda_1 = param(13);

% compute the non-local term
bb_nl = dot(llp,dmm)*phi*ones(nno,1);

% compute the reaction term
rt_1 = -cc+lambda_1*tau-llk.*(A*(1-tau)./(1+pp.^2).*...
       df_fun((1-tau).*kk,alpha_1,alpha_2,gamma)-delta_1-cc);
rt_1 = rt_1-(1-tau).*df_fun((1-tau).*kk,alpha_1,alpha_2,gamma).*theta.*bb_nl;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rt_2 = rtb2_fun(llk,llp,kk,pp,cc,tau,param)
%
% This functions constructs the reaction term of the equation for lambda_p(x,t)
%
% INPUT
%
% llk   := current lambda_k(x,t) distribution
% llp   := current lambda_p(x,t) distribution
% kk    := current state variable k(x,t) distribution
% pp    := current state variable p(x,t) distribution
% cc    := current control c(x,t) distribution
% tau   := current control tau(x,t) distribution
% param := parameters
%
% OUTPUT
%
% rt_2  := reaction term of the first equation evaluated in each mesh node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rt_2 = rtb2_fun(llk,llp,kk,pp,cc,tau,param)

% get the number of nodes
nno = size(kk,1);

% initialize the output
rt_2 = zeros(nno,1);

% extract useful parameters
delta_2 = param(4);
A = param(5);
alpha_1 = param(8);
alpha_2 = param(9);
gamma = param(10);
lambda_0 = param(12);

% compute the reaction term
rt_2 = lambda_0-llk.*(-2*A*pp./(1+pp.^2).^2.*...
       f_fun((1-tau).*kk,alpha_1,alpha_2,gamma))+delta_2*llp;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function beta_st = comp_beta_st(llk,llp,kk,pp,cc,tau,param,dmm)
%
% This functions compute beta_star(x,t) according to equation (1.14) 
%
% INPUT
%
% llk   := current lambda_k(x,t) distribution
% llp   := current lambda_p(x,t) distribution
% kk    := current state variable k(x,t) distribution
% pp    := current state variable p(x,t) distribution
% cc    := current control c(x,t) distribution
% tau   := current control tau(x,t) distribution
% param := parameters
% dmm   := diagonal mass matrix, obtained using mass lumping
%
% OUTPUT
%
% beta_st  := beta_star(x,t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta_st = comp_beta_st(llk,llp,kk,pp,cc,tau,param,dmm)

% get the number of nodes
nno = size(kk,1);

% initialize the output
beta_st = zeros(nno,1);

% extract useful parameters
delta_1 = param(3);
A = param(5);
theta = param(7);
alpha_1 = param(8);
alpha_2 = param(9);
gamma = param(10);
phi = param(11);
lambda_1 = param(13);

% compute the non-local term
bb_nl = dot(llp,dmm)*phi*ones(nno,1);

% compute the reaction term
beta_st = lambda_1+df_fun((1-tau).*kk,alpha_1,alpha_2,gamma).*...
          (A*llk./(1+pp.^2)+theta.*bb_nl);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function I = comp_I(llk,llp,kk,pp,cc,tau,param,dmm)
%
% This function computes the cost functional I 
%
% INPUT
%
% llk   := current lambda_k(x,t) distribution
% llp   := current lambda_p(x,t) distribution
% kk    := current state variable k(x,t) distribution
% pp    := current state variable p(x,t) distribution
% cc    := current control c(x,t) distribution
% tau   := current control tau(x,t) distribution
% param := parameters
% dmm   := diagonal mass matrix, obtained using mass lumping
%
% OUTPUT
%
% I     := cost functional I 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = comp_I(mkk,mpp,mcc,mtau,param,tt,dmm)

nt = length(tt);

It = zeros(1,nt);

lambda_0 = param(12);
lambda_1 = param(13);

for i = 1:nt
   cc = mcc(:,i);
   tau = mtau(:,i);
   kk = mkk(:,i);
   pp = mpp(:,i);

   It_1 = dot(dmm,cc.*kk);
   It_2 = dot(dmm,pp);
   It_3 = dot(dmm,tau.*kk);

   It(i) = It_1-lambda_0*It_2-lambda_1*It_3;
end

I = trapz(tt,It);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [node,elem] = get_mesh_P1_1d(x_a,x_b,nx)
% 
% This functions contsructs the data structure of a 1D mesh 
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
% node := node matrix ((nx+1)*1 matrix), containing the x-coordinates 
%         of the discretizations nodes 
% elem := elem (connectivity) matrix (nx*2), containing for each sub-interval, 
%         the indices of the two vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA,cc] = get_mat_P1_1d(node,elem,h)
%
% This functions constructs the stiffness and mass matrices 
% for 1D P1 finite elements
%
% INPUT
%
% node := node matrix ((nx+1)*1 matrix), containing the x-coordinates 
%         of the discretizations nodes
% elem := elem (connectivity) matrix (nx*2), containing for each sub-interval, 
%         the indices of the two vertices
% h    := mesh size (length of the generic sub-interval) 
%
% OUTPUT
%
% AA := stiffness matrix
% cc := diagonal mass matrix, obtained using mass lumping 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [AA_loc,cc_loc] = get_mat_P1_1d_loc(h)
%
% This functions constructs the local stiffness and mass matrices 
% for 1D P1 finite elements
%
% INPUT
%
% h    := mesh size (length of the generic sub-interval)
%
% OUTPUT
%
% AA_loc := local stiffness matrix
% cc_loc := local diagonal mass matrix, obtained using mass lumping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
