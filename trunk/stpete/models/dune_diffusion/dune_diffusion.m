function [hh] = dune_diffusion(x,h_to,nu)
% test the adams bashforth model
% try the case of a bump on flat bed
% let Q = 1/h
% formulate the solution as:
% dh/dt = h^(-2)dh/dx + nu d2h/dx2
% nu is a diffusion coef.

% initialize flat bottom
NX = length(x);
L = max(x)-min(x);
dx = x(2)-x(1);

% initialize h with bump or hole
figure;
plot(x, h_to)
pause(.5)

% now, setup equations
% constants
% stability for advective part
c_to = vfac*(h_to.^(-2));  
cmax = max(c_to);
dt_c = dx/cmax;

% stability for diffusive part
%nu = 1e-4;
dt_nu = 0.5*(dx^2)/max(nu);

% strictest stability requirement
dt = min([dt_c dt_nu]);

% factor of m for certain stability
m=5;
dt = dt/m;

% Witham defines a Reynolds type number
%A = dx*sum(c_to);
%Re = A/(2*nu);

% propose to run it for some amount of time
T = 6*L/(2*cmax);
NT = T/dt;

% initialize calculation arrays (save n_save times)
n_save = 100;
nt_save = ceil(NT/n_save);
hh = zeros(n_save,NX);
F = zeros(1,NX);

% initialize boundary conditions on h
h_last(1,:) = h_to;
hh(1,:) = h_last;

% initialize boundary conditions on F
% F is dQ/dx
% Q is 1/h
F(1) = 0;
F(NX) = 0;
F_last = F;

% begin
oo2dx = 1/(2*dx); % "one over delta-x"
oodx2 = 1/(dx*dx);  % "one over delta-x^2"
dQdx = zeros(1,NX);
nud2hdx2 = zeros(1,NX);
j_save = 0;
[~,imax] = max(h_last);
for j = 1:NT
    % find Dhigh?
    
    % spatial derivatives:
    for i= 2:(NX-1)
        h = h_last;
        dQdx(i) =  -( (1/h(i+1)) - (1/h(i-1)) ) * oo2dx;
        nud2hdx2(i) =  nu(i) * (h(i+1) - 2*h(i) + h(i-1)) * oodx2;
        F(i) = nud2hdx2(i);
    end
    
    % initialize F_last along array if i = 1
    if (j==1)
        F_last = F;
    end
    
    % do the time step
    h_last = dt * (1.5*F - 0.5*F_last) + h;
    F_last = F;
    
    % save results
    if ( (j-nt_save) > (j_save*nt_save))
        j_save = j_save + 1;
        hh(j_save,:) = h_last;
        fprintf('j = %d\r', j); 
    end
end

