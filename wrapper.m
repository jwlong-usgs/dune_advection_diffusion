% function [x, z, Dhigh, Dhighx, Dlow, Dlowx, counto, countc]=wrapper(xbliteGRIDS100,xb,xbliteHydro,xbliteD_all,profiles,nu)
%% WRAPPER
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Variables
% these variables only need to be loaded once. if you want to run different
% profiles, I would suggest evaluating the code in sections.
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\XB_HYDRO.mat'
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\XB_GRIDS100.mat'
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\xblite_parameters\xb_parameters.mat'
load '\\igsafpesvs002\StPetersburg-G_Shared\NACCH\Model\Data\Sandy_2012\XB_D_all.mat'
%% Convert Variables
% convert to model format
clearvars -except xbliteGRIDS100 xb profiles nu xbliteHydro xbliteD_all

profiles=2884;
% model coefficients
nuval=0.005;
vfac=10; % single value, may want to change

slope=0.52; % critical slope before avalanching
% hydro and features

t=3600:3600:180*3600;
twl=xbliteHydro.twlts(:,profiles)';
R2=xbliteHydro.R2ts(:,profiles)';
setup=xbliteHydro.setupts(:,profiles)';
S=xbliteHydro.Sts(:,profiles)';
surge=xbliteHydro.wlts(:,profiles)';
T=xbliteHydro.Tts(:,profiles)';
Ho=xbliteHydro.Hts(:,profiles)';
Bo=-xb.dlowslope(profiles,1); % foreshore slope: -xbliteD_all(profiles,13);

%% Models


Dhigh(:,1)=xb.preDhigh(profiles,1);
Dhighx(:,1)=xb.xpreDhigh(profiles,1);
Dlow(:,1)=xb.predlowf(profiles,1);
Dlowx(:,1)=xb.xpredlowf(profiles,1);
dVResidual=zeros(length(profiles),length(t));

for i=1:length(profiles)
    x(i).data(:,1)=(0:.1:((length(xbliteGRIDS100.pre.cXg(profiles).data)-1)/10))';
    z(i).data(:,1)=(xbliteGRIDS100.pre.cZi(profiles).data)';
    z(i).data(:,1)=flipud(z(i).data(:,1));
    nu=nuval.*ones(size(z(i).data(:,1))); % diffusion coefficient matrix size(z)
    counto(i)=0;
    countc(i)=0;
    zNewl=nan(length(t),length(z(1).data(:,1))); % test difference between coll/diff coll
%     Dlows=nan(length(t),length(z(1).data(:,1))); % change to structure if run on multiple profiles
    for j=1:length(t)
            % INUNDATION
        if surge(i,j)+setup(i,j)>=Dhigh(i,j) || isnan(Dhigh(i,j))
            z(i).data(:,j+1)=nan(size(z(i).data(:,1)));
            Dhigh(i,j+1)=NaN; % consider filling to 180, if you can bypass for loop...
            Dhighx(i,j+1)=NaN;
            Dlow(i,j+1)=NaN;
            Dlowx(i,j+1)=NaN;
            continue
        end
        
            % OVERWASH
        if twl(i,j)>=Dhigh(i,j) && surge(i,j)+setup(i,j)<Dhigh(i,j)
            [zNew] = dune_advection_diffusion(x(i).data(:,1),z(i).data(:,j),vfac,nu);
            [Dlowx(i,j+1), Dlow(i,j+1), Dhighx(i,j+1), Dhigh(i,j+1)]=find_dlow_dhigh(x(i).data(:,1),zNew);
            z(i).data(:,j+1)=zNew;
            counto(i)=counto(i)+1;
            
            
            % COLLISION
        elseif twl(i,j)>Dlow(i,j) && twl(i,j)<Dhigh(i,j)
            % if we want to run multiple profiles, will have to change
            % Dlows to a structure.
            [zNewl(j,:),dVResidual(i,j+1),Dlows(countc+1,:)] = LEH04_notime(x(i).data(:,1),z(i).data(:,j),...
                Dlow(i,j),Dlowx(i,j),3600,surge(i,j),T(i,j),Bo(i,1),R2(i,j),setup(i,j),S(i,j),dVResidual(i,j)); %add back Cs later
            
            [~,~,~,imin]=extreme(zNewl(j,:));
            % if imin is only one value, just run the entire profile
            if length(imin)<2
                imin(2)=length(x(i).data(:,1));
            end
            gridx=sort(imin);
            gridrx=gridx(2); % make sure this is pulling the second low
            [zNew] = dune_diffusion(x(i).data(:,1),zNewl(j,:),nu,vfac,gridrx); % slopec,
            [Dlowx(i,j+1), Dlow(i,j+1), Dhighx(i,j+1), Dhigh(i,j+1)]=find_dlow_dhigh(x(i).data(:,1),zNew,Dlows(1,:)');
            z(i).data(:,j+1)=zNew(:,1);
            countc(i)=countc(i)+1
            
            
            
            % SWASH
        else
            z(i).data(:,j+1)=z(i).data(:,j);
            Dhigh(i,j+1)=Dhigh(i,j);
            Dhighx(i,j+1)=Dhighx(i,j);
            Dlow(i,j+1)=Dlow(i,j);
            Dlowx(i,j+1)=Dlowx(i,j);
        end
    end
end

%% plot collision

figure;
subplot(2,1,1)
plot(Dlows(1,:))
hold on
plot(zNewl')
title(['nu = ' num2str(nuval)])

subplot(2,1,2)
plot(z(1).data)
hold on;
plot(fliplr(xbliteGRIDS100.pre.cZi(profiles).data),'--','LineWidth',2)
plot(fliplr(xbliteGRIDS100.post.cZi(profiles).data),'--','LineWidth',2)
plot(Dlows(1,:))

%% plot hydro
profile=1;
figure;
hold on; plot(surge(profile,:)+R2(profile,:))
plot(Dlow(profile,:))
plot(Dhigh(profile,:))
box on
grid on
title 'hydro vs time dependent dune'
legend('TWL','D_l_o_w','D_h_i_g_h')