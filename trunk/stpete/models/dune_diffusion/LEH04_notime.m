function [zNew,dVResidualn,Dlows] = LEH04_notime(x, z, Dlow, Dlowx, dt, surge, T, Bo, R2, setup, S, dVResidual,Cs)
% main program for revised LEH04 model with Palmsten and Holman 11
% updates.
% LEH04 - Larson, Erikson, Hanson (2004). An analyitical model to predict
% dune erosion due to wave impact. Coastal Eng. 51, 675-696.
% PH11a - Palmsten, Holman (submitted). Laboratory investigation of dune
% erosion using stereo video. Submitted to Coastal Eng.
%
%%  INPUTS:
%      x, z = cross-shore position and elevation, respectively
%      Dlowx, Dlow = cross-shore position and elevation of the Dune toe, respectively
%      t = time
%      surge = time series of storm surge
%       T = time series of peak wave period
%      B0 = whatever slope you want for the trajectory of the dune toe
%      R2 = 2% exceedence of wave runup
%      etabar = wave setup
%      sigma_s = swash (incident and infragravity combined)
%      Cs = sediment transport coefficient in LEH and Palmsten and Holman
%
%  OUTPUTS:
%      zNew = matrix of dune profiles through time
%
%%  Jwlong = revised code from Meg Palmsten and Kristen Splinter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 13
    Cs = 2.5e-4;
end

%% define constants
% dt = 3600;  % MAKE SURE THIS IS SECONDS % JACI took out t(2)-t(1), since t doesn't show up later
%nsigma = 2;   %in definition of R2, R16.. For R2 = nsigma=2, R16 = nsigma=1;
Kd = 1.26; % coefficient to account for higher runup on dune

% we can probably eliminate this.
Btfac = 1;
Bt = Bo*Btfac;

%%  Some hydrodynamic stuff
zTotal = R2.*Kd + surge;

% allocate variables

zNew = nan(length(x));
% if length(Cs) == 1
%     Cs = Cs(ones(size(x)));
% end


% %trajectory that dune toe receeds.
zb(1,1) = Dlow;
[~, st1] = min(((x)-Dlowx).^2);  %find grid point where initial dune toe is
% keyboard
zbT = [nan(st1-1,1); Bt.*(x(st1:end)-x(st1)) + zb(1)];

%% main program
% for tt=1:length(surge) % JACI just run one time, surge loop on outside
%     if tt==1
st = st1;
%     else
%
%     end
xToe = x(st);
V = sum(abs(diff(x(1:2))).*(z(st:end)));    %measured in ref to z=0
clear Vc
Vc = cumsum(abs(diff(x(1:2))).*(z(st:end)-zbT(st:end)));  %cumulative volume above the dune trajectory
Vc = Vc - Vc(1);

p =1-cdf('norm',zb,setup+surge,S);
Nc = p.*(dt./T);

%     if tt>1
dV = 4.*Cs.*(max(zTotal-zb,0)).^2.*Nc;
dVT = dV - dVResidual;
%     else
%         dVT = 4.*Cs.*(max(zTotal-zb,0)).^2.*Nc;
%     end
%
%     if dVT<0
%         ii=1;
%     else

[~, ii] = (min((Vc-dVT).^2)); %find grid point where dune toe is
%     end

%     dx = x(st+ii-1)-xToe;
dVResidualn = Vc(ii)-dVT;
%     zb(tt+1) = Bt.*dx + zb;  %trajectory that dune toe receeds.
st = st+ii-1;
zNew = [zbT(st).*ones(st,1); z(st+1:end)];
Dlows = Bt.*(x-x(st1)) + zb(1);



