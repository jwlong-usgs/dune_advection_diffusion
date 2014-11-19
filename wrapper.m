% function [x, z, Dhigh, Dhighx, Dlow, Dlowx]=wrapper(xbliteGRIDS100,xb,profiles,nu)
%% WRAPPER
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Variables
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\XB_HYDRO.mat'
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\XB_GRIDS100.mat'
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\xblite_parameters\xb_parameters.mat'

%% Convert Variables
% convert to model format

% model coefficients
nuval=0.05;
vfac=10; % single value, may want to change


% hydro and features
profiles=2884;
t=3600:3600:180*3600;
twl=xbliteHydro.twlts(:,profiles)';
R2=xbliteHydro.R2ts(:,profiles)';
setup=xbliteHydro.setupts(:,profiles)';
S=xbliteHydro.Sts(:,profiles)';
surge=xbliteHydro.wlts(:,profiles)';
T=xbliteHydro.Tts(:,profiles)';
Ho=xbliteHydro.Hts(:,profiles)';
Bo=-xb.dlowslope(profiles,1);

%% Models


Dhigh(:,1)=xb.preDhigh(profiles,1);
Dhighx(:,1)=xb.xpreDhigh(profiles,1);
Dlow(:,1)=xb.predlowf(profiles,1);
Dlowx(:,1)=xb.xpredlowf(profiles,1);
dVResidual=zeros(length(profiles),length(t));

for i=1:length(profiles)
    x(i).data(1,:)=0:.1:((length(xbliteGRIDS100.pre.cXg(profiles).data)-1)/10);
    z(i).data(1,:)=xbliteGRIDS100.pre.cZi(profiles).data;
    z(i).data(1,:)=fliplr(z(i).data(1,:));
    nu=nuval.*ones(size(z(i).data(1,:))); % diffusion coefficient matrix size(z)
    for j=1:length(t)
            % INUNDATION
        if surge(i,j)+setup(i,j)>=Dhigh(i,j) || isnan(Dhigh(i,j))
            z(i).data(j,:)=nan(size(z(i).data(1,:)));
            Dhigh(i,j+1)=NaN; % consider filling to 180, if you can bypass for loop...
            Dhighx(i,j+1)=NaN;
            Dlow(i,j+1)=NaN;
            Dlowx(i,j+1)=NaN;
            continue
        end
        
            % OVERWASH
        if twl(i,j)>=Dhigh(i,j) && surge(i,j)+setup(i,j)<Dhigh(i,j)
            [zNew] = dune_advection_diffusion(x(i).data(1,:),z(i).data(j,:),vfac,nu);
            [Dlowx(i,j+1), Dlow(i,j+1), Dhighx(i,j+1), Dhigh(i,j+1)]=find_dlow_dhigh(x(i).data(1,:),zNew);
            z(i).data(j+1,:)=zNew;
            
            % COLLISION
        elseif twl(i,j)>Dlow(i,j) && twl(i,j)<Dhigh(i,j)
            [zNew,dVResidual(i,j+1)] = LEH04MainProgram_v2(x(i).data(1,:),z(i).data(j,:),...
                Dlow(i,j),Dlowx(i,j),surge(i,j),T(i,j),Bo(i,1),R2(i,j),setup(i,j),S(i,j),dVResidual(i,j)); %add back Cs later
            
            [zNew] = dune_diffusion(x(i).data(1,:),zNew,nu,vfac);
            [Dlowx(i,j+1), Dlow(i,j+1), Dhighx(i,j+1), Dhigh(i,j+1)]=find_dlow_dhigh(x(i).data(1,:),zNew);
            z(i).data(j+1,:)=zNew;
            
            % SWASH
        else
            z(i).data(j+1,:)=z(i).data(j,:);
            Dhigh(i,j+1)=Dhigh(i,j);
            Dhighx(i,j+1)=Dhighx(i,j);
            Dlow(i,j+1)=Dlow(i,j);
            Dlowx(i,j+1)=Dlowx(i,j);
        end
    end
end