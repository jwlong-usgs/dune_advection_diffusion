
function [Dlowx, Dlow, Dhighx, Dhigh]=find_dlow_dhigh(x,z,Dlows)

%inputs:
% x,z profile cross-shore position and elevation

% outputs:
% Dlowx = cross-shore position of Dlow
% Dlow = elevation of Dlow
% Dhighx = cross-shore position of Dhigh
% Dhigh = elevation of Dhigh


%%%%%%%%%%%%% find Dhigh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xmax,imax,~, ~] = extreme(z);

if imax == 1
    xmax = xmax(2:end);
    imax = imax(2:end);
end

Dhighx = x(imax(1));
Dhigh = xmax(1);

%%%%%%%%%%%%% find Dlow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xb = x(x <= Dhighx);
zb = z(x <= Dhighx); 
Dlowsb = Dlows(x <= Dhighx);
% allocate
% dzdx = nan*ones(length(xb),1);
% 
% for jj= 2:length(xb)-1
%     
%     dzdx(jj) = (zb(jj-1)-zb(jj+1))/(xb(jj-1)-xb(jj+1)); % endpoint slope
%     
% end
% 
% ddz = diff(dzdx); %find where 2nd derivative is maximum
% [~,ind] = max(ddz);
% 
% Dlow = mean(zb(ind:ind+1));
% Dlowx = mean(xb(ind:ind+1));
% keyboard
[r,~] = find(zb>Dlowsb);
Dlowx = xb(r(1),1);
Dlow = zb(r(1),1);
