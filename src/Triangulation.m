function [X_Y_Z_Model,var_py,py]=Triangulation(unknown0,XY_L,XY_R,f,xp,yp,k1,k2,p1,p2)
% XY_L=Points_Im_1_Image_Cs;XY_R=Points_Im_2_Image_Cs;
n=size(XY_R,1);
Qy=eye(4*n);
%-----------------In-----------
r_L=((XY_L(:,1)-xp.*ones(n,1)).^2 + (XY_L(:,2)-yp.*ones(n,1)).^2).^0.5;
dx_L=( k1.*r_L.^2+k2.*r_L.^4).*(XY_L(:,1)-xp.*ones(n,1))+p2.*(r_L.^2+2.*(XY_L(:,1)-xp.*ones(n,1)).^2)+2*p1.*(XY_L(:,1)-xp.*ones(n,1)).*(XY_L(:,2)-yp.*ones(n,1));
dy_L=( k1.*r_L.^2+k2.*r_L.^4).*(XY_L(:,2)-yp.*ones(n,1))+p1.*(r_L.^2+2.*(XY_L(:,2)-yp.*ones(n,1)).^2)+2*p2.*(XY_L(:,1)-xp.*ones(n,1)).*(XY_L(:,2)-yp.*ones(n,1));
r_R=((XY_R(:,1)-xp.*ones(n,1)).^2 + (XY_R(:,2)-yp.*ones(n,1)).^2).^0.5;
dx_R=( k1.*r_R.^2+k2.*r_R.^4).*(XY_R(:,1)-xp.*ones(n,1))+p2.*(r_R.^2+2.*(XY_R(:,1)-xp.*ones(n,1)).^2)+2*p1.*(XY_R(:,1)-xp.*ones(n,1)).*(XY_R(:,2)-yp.*ones(n,1));
dy_R=( k1.*r_R.^2+k2.*r_R.^4).*(XY_R(:,2)-yp.*ones(n,1))+p1.*(r_R.^2+2.*(XY_R(:,2)-yp.*ones(n,1)).^2)+2*p2.*(XY_R(:,1)-xp.*ones(n,1)).*(XY_R(:,2)-yp.*ones(n,1));
XY_L(:,1)=XY_L(:,1)+dx_L;
XY_L(:,2)=XY_L(:,2)+dy_L;
XY_R(:,1)=XY_R(:,1)+dx_R;
XY_R(:,2)=XY_R(:,2)+dy_R;
%----------------In------------
omega2=unknown0(1);phi2=unknown0(2);kapa2=unknown0(3);by=unknown0(4);bz=unknown0(5);bx=unknown0(6);
% disp('     omega2                phi2                    kapa2                  by                 bz')

%Rotation matrix
R_kapa = [cos(kapa2) sin(kapa2) 0 ; -sin(kapa2) cos(kapa2) 0 ; 0 0 1];
R_phi = [cos(phi2) 0 -sin(phi2) ; 0 1 0 ; sin(phi2) 0 cos(phi2)];
R_omega = [1 0 0 ; 0 cos(omega2) sin(omega2) ; 0 -sin(omega2) cos(omega2)];
R = R_kapa*R_phi*R_omega;
%Triangulation-------
xxx=[XY_L(:,1)-xp XY_L(:,2)-yp -f*ones(n,1)];
x_1=xxx(:,1);
y_1=xxx(:,2);
z_1=xxx(:,3);
XXX=[XY_R(:,1)-xp XY_R(:,2)-yp -f*ones(n,1)]*R';
x_2=XXX(:,1);
y_2=XXX(:,2);
z_2=XXX(:,3);
landa=(bx.*z_2-bz.*x_2)./(x_1.*z_2-x_2.*z_1);
mu=(bx.*z_1-bz.*x_1)./(x_1.*z_2-x_2.*z_1);
%Model Coordinates
Xm=landa.*x_1;
Zm=landa.*z_1;
Ym=0.5*landa.*y_1+0.5*mu.*y_2+0.5*by;
py=landa.*y_1-mu.*y_2-by;
var_py=sqrt(mean(py.^2));
X_Y_Z_Model=[Xm,Ym,Zm];
end