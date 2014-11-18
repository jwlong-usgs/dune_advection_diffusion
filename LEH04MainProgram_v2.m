function [xM,zNew,dVT] = LEH04MainProgram_v2(xzGrid,xzFinal,t,WL,Ho,T,dLslope)
%% main program for revised LEH04 model with Palmsten and Holman 11
%% updates.
%% LEH04 - Larson, Erikson, Hanson (2004). An analyitical model to predict
%% dune erosion due to wave impact. Coastal Eng. 51, 675-696.
%% PH11a - Palmsten, Holman (submitted). Laboratory investigation of dune
%% erosion using stereo video. Submitted to Coastal Eng.

%addpath C:\Users\mpalmsten\Documents\OSUStuff\matlab\classes\wavesClass\project\matlab\
%close all
%clear all

%load data inital profile
%xzGrid = load('MTA67059xz.dat');
xTemp = xzGrid(:,1);
zTemp = xzGrid(:,2);
clear xzGrid;

% clip bathy deeper than 1 m
keep = find(zTemp>=-.25);
x = xTemp(keep);
z = zTemp(keep);
clear xTemp zTemp

% load final profile
%xzFinal = load('MTA67060xz.dat');
zTemp = xzFinal(:,2);
xTemp = xzFinal(:,1);
keep = find(zTemp>=-1);
zFinal = zTemp(keep);
xFinal = xTemp(keep);

clear zxFinal

% interpolate data to a uniform grid and start with origin at the offshore
% boundary 
xF = xTemp(keep);
clear zTemp xTemp
shift = xF(end)-x(end)
xM = [0:1:max(xF(end)-min(xF), x(end)-min(x))]; %put this on a uniform grid % [1262:.1:1350]; JACI
% xM = [0:1:max(xF(1)-min(xF), x(1)-min(x))];% changed to x(1), because the profiles are backwards before being put into the model
zTemp = z;
clear z;
z = interp1(max(x)-x,(zTemp),xM);
clear zTemp
zTemp=zFinal;
clear zFinal
zFinal = interp1(max(x)-xF,zTemp,xM);
xM = xM % -1262;JACI
%waves and water levels
% WaveData = load('waveConditionsGoldCoast.txt');
% WLData = load('waterLevelGoldCoast.txt');
% t1 = WLData(:,1);
% dt1 = diff(t1(1:2));
% WL = WLData(:,2);
% Hs = WaveData(:,1);
% T = WaveData(:,2);
% theta = WaveData(:,3);  %need to account for wave refraction, so need to get offshore wave height equivalent...
% dt2 = WaveData(1,6);
% t2 = [0:dt2:(length(Hs)-1).*dt2];
% if(dt1~=dt2)
%     dt = max(dt1,dt2);
%      t = [0:dt:t2(end)];
%     if (dt==dt1)
%         HsTemp = Hs;
%         TTemp = T;
%         clear Hs T
%         Hs = interp1(t2,HsTemp,t);
%         T = interp1(t2,TTemp,t);
%     else
%         WLTemp = WL;
%         clear WL
%         WL = interp1(t1,WLTemp,t);
%     end
% end
% %% SHOULD THIS BE Ho,rms?
% Ho = calcHoFromHs(Hs,T,20.*ones(length(T),1)); %note that Hill's model reverse shoals and doesn't account for refraction
Hrmso = Ho./sqrt(2);

dt = t(2)-t(1);  % MAKE SURE THIS IS SECONDS

%% attempts to find zb(0) and Btfac;
figure
plot(xM,z)
grid on
clear keep
%ch = ginput(2)
[val zS] = (min((z-1).^2)); %used to find swash/beach slope --> how does this define the swash zone?
[val dC] = max(z); %Dune crest


% this attempt at finding the dune toe does not appear to have worked
% keep = find(xM>xM(zS)-5 & xM<xM(dC));
% B1 = regress(z(keep)',[ones(size(z(keep)')) xM(keep)'])
% hold on, plot(xM(keep),([ones(size(z(keep)')) xM(keep)']*B1),'g')
% clear keep
% keep = find(xM<xM(zS)+5);
% B2 = regress(z(keep)',[ones(size(z(keep)')) xM(keep)']);
% hold on, plot(xM(keep),([ones(size(z(keep)')) xM(keep)']*B2),'g')
% cubicFit = 3.1e-005.*xM.^3 - 0.0054.*xM.^2 + 0.33.*xM - 7.2
% keep=find(xM>50 &xM<110);
% [val,ii]=max(cubicFit(keep) - z(keep))
% zb2 = z(keep(ii))
% hold on, plot(xM,zFinal,'r')
%z2 = (z(:)-[ones(size(z')) xM']*B2).^2
%z3 = (z(:)-[ones(size(z')) xM']*B1).^2
%zR = min(z2,z3)
%ch2=ginput(2)
%ch2 = [68.2258 2.1579;
%  133.3410 -0.8246];
%clear keep
%keep = find(xF>ch2(1)& xF<ch2(2));
%B3 = regress(zFinal(keep),[ones(size(zFinal(keep))) xF(keep)'])
%hold on, plot(xF(keep),([ones(size(zFinal(keep))) xF(keep)']*B3),'r')

%% constants
nsigma = 2; %in definition of R2, R16.. For R2 = nsigma=2, R16 = nsigma=1; 
g = 9.81;
Bo = dLslope; %3.5./(55);  %initial beach slope, between dune toe and z=0, Dec. 2008 fig.
Btfac = -1; %(2.183-.9471)./(112-85)./Bo;  %slope at which beta receeds. LEH04 = 1, PH11 = 0.54.... 
Bt = Bo*Btfac;
Bf = 2.183./(112-79);   %final beach slope
Bf = .5/30
zb(1) =3.03681683881661;   %based on DEc. 2008 survey, regression fit. See fig.
[val st1] = (min((z-zb(1)).^2)); %find grid point where initial dune toe is
st1=550; 
Cs = 1.7*10^(-4).*ones(size(Hrmso));  %in LEH04 model, 2.5*10^-4 is Birk data set max, LAB data 1.4 x 10^-3
Cs = 20*10^(-3).*ones(size(Hrmso));
%Cs = 20*10^(-3);    %best fit for ETA67 using Meg's model, but R2
%Cs = 500*10^(-3);   %still not enough to get good match using R16
Ac = 1.34*10^(-3);
bc = 3.19*10^(-4);
D50 = 0.000201; % from Gold Coast Shoreline Mangagement Plan
%Cs = Ac.*exp(-bc.*Hrmso./D50);
Kd = 1.26; % coefficient to account for higher runup on dune

zbT = [repmat(NaN,1,st1-1) Bt.*(xM(st1:end)-xM(st1)) + zb(1)];  %trajectory that dune toe receeds.
hold on, plot(xM,zbT,'k-')
Lo = g.*T.^2./(2.*pi);
clear ii tt
%% main program
 for tt=1:length(WL)
%[val st] = (min((z-zb(tt)).^2)); %find grid point where dune toe is
if tt==1
    st = st1;
else
    st = ii;
end
xToe(tt) = xM(st);
V(tt) = sum(abs(diff(xM(1:2))).*(z(st:end)));    %measured in ref to z=0
clear Vc
 Vc = cumsum(abs(diff(xM(1:2))).*(z(st:end)-zbT(st:end)));  %cumulative volume above the dune trajectory
 Vc = Vc - Vc(1);
 Beta(tt) = Bo;   %initial dummy guess.
 etabar(tt) = 0.35.*Beta(tt).*sqrt(Ho(tt).*Lo(tt));
 sigma_s(tt) = sqrt(Ho(tt).*Lo(tt).*(0.563.*(Beta(tt).^2)+0.0004))./2.*nsigma./2;
 zR(tt) = 1.1.*(etabar(tt)+ sigma_s(tt));
 sigma_s2(tt) = sqrt(Ho(tt).*Lo(tt).*(0.563.*(Beta(tt).^2)+0.0004))./2;
 zR2(tt) = 1.1.*(etabar(tt)+ sigma_s2(tt));
 zRLEH(tt) = 0.158.*sqrt(Ho(tt)./1.416.*Lo(tt)); 
 %zTotal(tt) = zRLEH(tt) + WL(tt);
 zTotal(tt) = zR2(tt).*Kd + WL(tt);
% 
%y = pdf('norm',[-2:0.01:2],etabar(tt),sigma_s(tt));

p(tt) =1-cdf('norm',zb(tt),etabar(tt)+WL(tt),sigma_s(tt));
 Nc(tt) = p(tt).*(dt./T(tt));
 if tt>1
        dV(tt) = 4.*Cs(tt).*(max(zTotal(tt)-zb(tt),0)).^2.*Nc(tt);
        dVT(tt) = dV(tt) - dVResidual(tt-1);
 else
     dVT(tt) = 4.*Cs(tt).*(max(zTotal(tt)-zb(tt),0)).^2.*Nc(tt);
 end
 %these cause a bit of a problem when the dune is quite large compared to
 %the runup... We'd still expect undercutting and failure

 %try estimating dx based on eroded volume
%Vc = cumsum(abs(diff(xM(1:2))).*(max(min(zTotal(tt),z(st:end))-zbT(st:end),0)));  %cumulative volume above the dune trajectory
%Vc = Vc - Vc(1);
xstep = 0;  dVtemp = 0; if tt==1; zcheck = z; dV(1) = 0; else zcheck = zNew(tt-1,:); end
imax = [];
if dVT(tt)<0
    ii=1;
else
    while dV(tt)>dVtemp
    zi = interp1(xM+xstep,z,xM);
    ind = find(zcheck>zb(tt));
    [dVtemp,imax] = max(cumtrapz(xM(zcheck>zb(tt)),zcheck(zcheck>zb(tt))-zi(zcheck>zb(tt))));
    imax = imax+ind(1)-2;
    xstep = xstep+.1; 
    end
[~, ii] = (min((xM-(xToe(tt)+xstep-.1)).^2)); %find grid point where dune toe is
%ii=ii-st1+1;
end

if isempty(imax); imax=st; zi = z; end
dx(tt) = xstep-.1; %xM(st+ii-1)-xToe(tt);
xToe(tt+1) = xToe(tt)+dx(tt);
dVResidual(tt) = 0; %Vc(ii)-dVT(tt);
zb(tt+1) = Bt.*dx(tt) + zb(tt);  %trajectory that dune toe receeds.
zforeslope = zb(tt+1)/(ii/10);
zfore = xM(1:ii)*zforeslope;
zNew(tt,:) = [zfore zi(ii+1:imax+1) z(imax+2:end)];
figure(111)
plot(xM,z,'-b')
hold on; pause
plot(xM,zNew(tt,:),'--r')
end
 
 figure(111); plot(xM,zFinal,'--k','linewidth',2)
%keyboard
xToe(tt+1) = xToe(tt)+dx(tt);

 st = st+ii-1;
 %zNew(tt+1,:) = [z(1:st1) zbT(st1+1:st) z(st+1:end)];
 t = [t t(end)+dt];
 figure
subplot(211),plot(t./60./60,zb,'k','linewidth',3)
hold on, ch=line([t(1) t(end)]./3600, [zb(1) zb(1)]);
set(ch,'linewidth',3)
hold on, ch=line([t(1) t(end)]./3600, [2.183 2.183]);
set(ch,'linewidth',3,'color','r')
grid on
legend('z_b(t)','z_b(0)', 'z_b(final) obs')
axis tight
ylim([0.5 2.5])
xlabel('t (hrs)')
ylabel('z (m)')
title('ETA 67')
subplot(212),plot(t./60./60,xToe,'k','linewidth',3)
hold on, ch=line([t(1) t(end)]./3600, [xToe(1) xToe(1)]);
set(ch,'linewidth',3)
hold on, ch=line([t(1) t(end)]./3600, [112 112]);
set(ch,'linewidth',3,'color','r')
grid on
legend('x_b(t)','x_b(0)', 'x_b(final) obs')
axis tight
ylim([80 120])
xlabel('t (hrs)')
ylabel('x (m)')

figure
plot(xM,z,'k','linewidth',3)
hold on, plot(xM,zFinal,'k--','linewidth',3)
hold on, plot(xM,zNew(end,:),'r')
hold on, ch=line([0 65],[0 0]);
set(ch,'linewidth',3)
legend('Dec 2008 profile', 'June 2009 profile','z(t)','SWL')
hold on, plot(xM,zNew(1:12:end,:),'r')
axis tight
ylim([-2 10])
xlabel('x (m)')
ylabel('z (m)')
title('ETA 67')

figure, plot(t(1:end-1)./3600,zTotal,'linewidth',3)
hold on, plot(t(1:end-1)./3600,zR+WL,'r','linewidth',3)
hold on, ch=line([t(1) t(end-1)]./3600,[zb(1) zb(1)])
set(ch,'lineStyle','--', 'color','k','linewidth',3)
grid on
axis tight
ylim([0 7])
legend('R_{LEH} + tide + surge','R_{2} + tide + surge','z_b(0)')
ylabel('z (m)')
xlabel('t (hrs)')

figure
plot(t(1:end-1)./3600,Cs,'linewidth',3)
grid on
axis tight
ylabel('C_s')
xlabel('t (hrs)')
