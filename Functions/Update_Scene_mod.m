function [d_n,theta_n,vf] = Update_Scene_mod(d,theta,v,ToS)
d_n = sqrt(d.^2 + (v.^2)*(ToS.^2)-2.*d.*v.*cosd(theta)*ToS);
theta_n = sign(theta).*abs(atand((d.*sind(theta))./(d.*cosd(theta)-v*ToS)));
vf = theta_n<90 | theta_n>90;