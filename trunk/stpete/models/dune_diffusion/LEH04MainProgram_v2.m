function [zNew] = LEH04MainProgram_v2(x, z, Dlow, Dlowx, t, surge, T, Bo, R2, etabar, sigma_s, Cs)
%% main program for revised LEH04 model with Palmsten and Holman 11
%% updates.
%% LEH04 - Larson, Erikson, Hanson (2004). An analyitical model to predict
%% dune erosion due to wave impact. Coastal Eng. 51, 675-696.
%% PH11a - Palmsten, Holman (submitted). Laboratory investigation of dune
%% erosion using stereo video. Submitted to Coastal Eng.
%%
%%  INPUTS:
%%      x, z = cross-shore position and elevation, respectively
%%      Dlowx, Dlow = cross-shore position and elevation of the Dune toe, respectively
%%      t = time
%%      surge = time series of storm surge
%       T = time series of peak wave period
%%      B0 = whatever slope you want for the trajectory of the dune toe
%%      R2 = 2% exceedence of wave runup
%%      etabar = wave setup
%%      sigma_s = swash (incident and infragravity combined)
%%      Cs = sediment transport coefficient in LEH and Palmsten and Holman
%%
%%  OUTPUTS:
%%      zNew = matrix of dune profiles through time
%%
%%  Jwlong = revised code from Meg Palmsten and Kristen Splinter

if nargin < 12
    Cs = 2.5e-4;
end

%% define constants
dt = t(2)-t(1);  % MAKE SURE THIS IS SECONDS
%nsigma = 2;   %in definition of R2, R16.. For R2 = nsigma=2, R16 = nsigma=1;
Kd = 1.26; % coefficient to account for higher runup on dune

% we can probably eliminate this.
Btfac = 1;
Bt = Bo*Btfac;

%%  Some hydrodynamic stuff
zTotal = R2.*Kd + surge;

% allocate variables
xToe = nan(size(surge));
V = nan(size(surge));
dV = nan(size(surge));
dVT = nan(size(surge));
dx = nan(size(surge));
dVResidual = nan(size(surge));
p = nan(size(surge));
Nc = nan(size(surge));
zb = nan(size(surge));
zNew = nan(length(surge),length(x));
if length(Cs) == 1
    Cs = Cs(ones(size(x)));
end


% %trajectory that dune toe receeds.
zb(1) = Dlow;
[~, st1] = (min((x-Dlowx).^2));  %find grid point where initial dune toe is
zbT = [nan(st1-1,1); Bt.*(x(st1:end)-x(st1)) + zb(1)];

%% main program
for tt=1:length(surge)
    if tt==1
        st = st1;
    else
        st = st+ii-1;
    end
    xToe(tt) = x(st);
    V(tt) = sum(abs(diff(x(1:2))).*(z(st:end)));    %measured in ref to z=0
    clear Vc
    Vc = cumsum(abs(diff(x(1:2))).*(z(st:end)-zbT(st:end)));  %cumulative volume above the dune trajectory
    Vc = Vc - Vc(1);
    
    p(tt) =1-cdf('norm',zb(tt),etabar(tt)+surge(tt),sigma_s(tt));
    Nc(tt) = p(tt).*(dt./T(tt));
    
    if tt>1
        dV(tt) = 4.*Cs(tt).*(max(zTotal(tt)-zb(tt),0)).^2.*Nc(tt);
        dVT(tt) = dV(tt) - dVResidual(tt-1);
    else
        dVT(tt) = 4.*Cs(tt).*(max(zTotal(tt)-zb(tt),0)).^2.*Nc(tt);
    end
    
    if dVT(tt)<0
        ii=1;
    else
        [~, ii] = (min((Vc-dVT(tt)).^2)); %find grid point where dune toe is
    end
    
    dx(tt) = x(st+ii-1)-xToe(tt);
    dVResidual(tt) = Vc(ii)-dVT(tt);
    zb(tt+1) = Bt.*dx(tt) + zb(tt);  %trajectory that dune toe receeds.
    
    zNew(tt,:) = [z(1:st1); zbT(st1+1:st); z(st+1:end)];
end




