function [ inival ] = initialval3dconf( xyz,XYZ )
%
%

confptr = twoDconformal(xyz(:,2:3),XYZ(:,2:3));

%initial val of Scale
L = sqrt((confptr(1)^2)+(confptr(2)^2));

%initial val of Omega and Phi rot.
Omg= 0;
Phi= 0;

%initial val of Kappa rot.
Kap= atan2(confptr(2),confptr(1));

MW = [1 0 0;0 cos(Omg) -sin(Omg);0 sin(Omg) cos(Omg)];
MP = [cos(Phi) 0 sin(Phi);0 1 0;-sin(Phi) 0 cos(Phi)];
MK = [cos(Kap) -sin(Kap) 0;sin(Kap) cos(Kap) 0;0 0 1];
M1 = MK*MP*MW;

n = size(xyz,1);
for j=1:n
    
    T = XYZ(j,2:4)' - (M1'/L)*xyz(j,2:4)';
    Translation(:,j) = T;
    
end
Translation = mean(Translation,2);

inival = [Omg Phi Kap Translation' L]';

end

