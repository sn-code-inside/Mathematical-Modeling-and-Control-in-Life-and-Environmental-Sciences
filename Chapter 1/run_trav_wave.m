%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAM RUN_TRAV_WAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program solves the controlled traveling wave system,
% equation (1.23),
% with the numerical methods described in Section 1.5.3,
% on the space domain [x_a,x_b]
% and the time interval [0,Tend].
% With the current set of parameters, it reproduces in particular
% the results of experiment 2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a = -60;
x_b = 2;
nx = 2000;
h = (x_b-x_a)/nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tend = 200;
dt = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1;     % diffusion coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MESH PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1 = nx+1;
nno = nx1;
nel = nx;
dofs = 2*nno;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE NODE AND ELEM MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[node,elem] = get_mesh_P1_1d(x_a,x_b,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE STIFFNESS AND LUMPED MASS MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA,cc] = get_mat_P1_1d(node,elem,h);
AA = d*AA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE ITERATION MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccs = spdiags(cc,0,nno,nno);
mZero = sparse(nno,nno);
mat = [1/dt*ccs+AA mZero;mZero 1/dt*ccs];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET INITIAL CONDITIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu_1 = zeros(nno,1);
uu_2 = zeros(nno,1);
for i = 1:nno
   x = node(i);
   if(x<=-1)
     uu_1(i) = 0;
     uu_2(i) = 0;
   elseif(x>-1&&x<=1)
     uu_1(i) = 1/2*(1+x);
     uu_2(i) = 1/2*(1+x);
   else
     uu_1(i) = 1;
     uu_2(i) = 1;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET CONTROL BETA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = zeros(nno,1);
for i=1:nno
   x = node(i);
   if(x<=-20)
     beta(i) = 2.1;
   else
     beta(i) = 1;
   end
end 
%beta = ones(nno,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET TIME PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = [0:dt:Tend];
nt = length(tt);
t = 0;
kk = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE SOLUTION MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muu = zeros(2*nno,ceil(nt/100));
muu(1:nno,1) = uu_1;
muu(nno+1:end,1) = uu_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START TIME LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:nt-1
  t = t+dt;

  % compute reaction terms
  bbr1 = rt1_fun(node,cc,uu_2);
  bb1 = cc.*(1/dt*uu_1-uu_1+bbr1);
  gg2 = g_fun(uu_1);
  bb2 = cc.*(1/dt*uu_2-beta.*uu_2+gg2);
  bb = [bb1;bb2];

  % solve the linear system
  uu = mat\bb;

  if(mod(k,100)==1)
     kk = kk+1;
     muu(:,kk) = uu;
     disp(['t = ',num2str(t)])
  end

  % update
  uu_1 =uu(1:nno);
  uu_2 =uu(nno+1:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_plot(muu,nno)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function gg = g_fun(uu_1)
%
% This functions evaluates g(u_1(x,t))
%
% INPUT
%
% uu_1  := current state variable u_1(x,t) distribution
%
% OUTPUT
%
% gg    := g_(u_1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gg = g_fun(uu_1)

gg=2*uu_1./(1+uu_1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bbr1 = rt1_fun(node,cc,uu_2);
%
% This functions constructs the reaction term of the equation for u_1(x,t)
%
% INPUT
%
% node  := node matrix ((nx+1)*1 matrix),
%          containing the x-coordinates of the discretizations nodes
% cc    := diagonal mass matrix, obtained using mass lumping
% uu_2  := current state variable u_2(x,t) distribution
%
% OUTPUT
%
% rt_1  := reaction term of the equation for u_1(x,t) evaluated in each mesh node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bbr1 = rt1_fun(node,cc,uu_2)

nno = size(uu_2,1);

bbr1 = zeros(nno,1);

sigma = 2;

for i = 1:nno
   x1 = node(i);
   kk = zeros(nno,1);
   for j=1:nno
      x2 = node(j);
      kk(j) = 1/sqrt(2*pi)*exp(-(x1-x2)^2/sigma);
   end
   bbr1(i) = dot(kk,cc.*uu_2);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [node,elem] = get_mesh_P1_1d(x_a,x_b,nx)
% 
% This function constructs the data structure of a 1D mesh 
% for P1 finite elements
%
% INPUT
% 
% x_a  := left vertex of the 1D interval domain
% x_b  := right vertex of the 1D interval domain
% nx   := number of 1D finite elements (sub-intervals)
%
% OUTPUT
%
% node := node matrix ((nx+1)*1 matrix), 
%         containing the x-coordinates of the discretizations nodes 
% elem := elem (connectivity) matrix (nx*2), containing for each sub-interval, 
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
% elem := elem (connectivity) matrix (nx*2), containing for each sub-interval, 
%         the indices of the two vertices
% h    := mesh size (length of the generic sub-interval) 
%
% OUTPUT
%
% AA   := stiffness matrix
% cc   := diagonal mass matrix, obtained using mass lumping 
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
% h      := mesh size (length of the generic sub-interval)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_plot(muu,nno)
%
% Plot spatio-temporal distribution of u_1(x,t) 
%
% INPUT
%
% muu    := spatio-temporal distributions of the state variables 
% nno    := number of nodes 
%
% OUTPUT
%
% AA_loc := local stiffness matrix
% cc_loc := local diagonal mass matrix, obtained using mass lumping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_plot(muu,nno)

nt = 201;
nx = 201;

tspan = linspace(0,200,nt);
xspan = linspace(-60,2,nx);

surf(xspan,tspan,muu(1:10:nno,:)')
axis([-60 2 0 200 0 1])
set(gca,'XTick',[-60:20:0],'YTick',[0:50:200],'ZTick',[0:0.5:1],'Fontsize',30)
xlabel('x','Fontsize',30)
ylabel('t','Fontsize',30)

view(-30,10)

end
