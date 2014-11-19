% setup for LEH
%
%
%
%

profile = 1;
endtime = 80; % time in hours, stop before overwash?
%% convert variables
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\xbtest.mat'
xzGrid(:,1)=0:.1:((length(xbtest.pre.cZi(profile).data)-1)/10);
%xzGrid = flipud(xzGrid);
%xzGrid(:,1)=max(xzGrid(:,1))-xzGrid(:,1);
xzGrid(:,2)=xbtest.pre.cZi(profile).data;
xzGrid(:,2) = flipud(xzGrid(:,2));


xzFinal(:,1)=xzGrid(:,1);
xzFinal(:,2)=xbtest.post.cZi(profile).data;
xzFinal(:,2)=flipud(xzFinal(:,2));

t=3600:3600:endtime*3600;
WL=xbtest.twl(1:endtime,profile)';
Ho=xbtest.H(1:endtime,profile)';
T=xbtest.T(1:endtime,profile)';
R2=xbtest.R2(:,profile)';
etabar=xbtest.setup(:,profile)';
sigma_s=xbtest.S(:,profile)';

dLslope=-xbtest.dlowslope(profile,1); %xbtest.dlowslope(profile,1);
bslope=xbtest.prebeachslope(profile,1);
%dLslope = -bslope; 
%% run model

%[xM,zNew,dVT] = LEH04MainProgram_v2(xzGrid,xzFinal,t,WL,Ho,T,dLslope,bslope,R2,etabar,sigma_s);

% for the new version
Dlowx = xbtest.Dlowfront(1,1);
Dlow = xbtest.Dlowfront(1,2);

[zNew] = LEH04MainProgram_v2(xzGrid(:,1), xzGrid(:,2), Dlow, Dlowx, t, WL, T, dLslope, R2, etabar, sigma_s);