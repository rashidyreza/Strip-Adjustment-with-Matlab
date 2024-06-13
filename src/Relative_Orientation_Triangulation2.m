function [X_Y_Z_Model,var_py,py,pixel_error,parameters,unknown0]=Relative_Orientation_Triangulation(XY_L,XY_R,f,xp,yp,Pixel_Size,k1,k2,p1,p2)



% xp=0;yp=0;
% [n1,n2]=uigetfile('*.txt','xyL');%xy photo grametric aks chap
% XY_L=strcat(n2,n1);
% XY_L=dlmread(XY_L);
% XY_L=0.001.*XY_L;
% 
% [n1,n2]=uigetfile('*.txt','xyR');%xy photo grametric aks rast
% XY_R=strcat(n2,n1);
% XY_R=dlmread(XY_R);
% XY_R=0.001.*XY_R;
% n=size(XY_R,1);
% Qy=eye(4*n);
% f = .152844;

%focal lens







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
syms omega2 phi2 kapa2 by bz %Unknown

%Rotation matrix
R_kapa = [cos(kapa2) sin(kapa2) 0 ; -sin(kapa2) cos(kapa2) 0 ; 0 0 1];
R_phi = [cos(phi2) 0 -sin(phi2) ; 0 1 0 ; sin(phi2) 0 cos(phi2)];
R_omega = [1 0 0 ; 0 cos(omega2) sin(omega2) ; 0 -sin(omega2) cos(omega2)];
R = R_kapa*R_phi*R_omega;

syms xL yL xR yR %Obz

bx = 0.01;
% Coplanarity Equation

F0=det( [bx by bz ; xL-xp yL-yp -f ;[xR-xp yR-yp -f]*R'] );

L = [xL yL xR yR];
L0=[XY_L(:,1),XY_L(:,2),XY_R(:,1),XY_R(:,2)]';
L0=reshape(L0,[],1);
L00=L0;
unknown = [omega2 phi2 kapa2 by bz];
B0 = jacobian(F0,L);
A0 = jacobian(F0,unknown);
unknown0 = [0 0 0 0 0];
e0=2;
e1=3;
dx=1;
B = zeros(n,4*n);
A = zeros(n,5);
W = zeros(n,1);
digits(100);
while (((norm(dx)< 1e-10)&&(norm(e1-e0))<1e-1))==0
    
    for i=1:n
    c = 4*i-3:4*i;
    a=A0;
    a = subs(a,L,L0(c)');
    a = subs(a,unknown,unknown0);
    a = eval(a);
    A(i, :) = a;
    

        b = B0;
        b = subs(b, unknown, unknown0);
        b = subs(b,L,L0(c)');
        b = eval(b);
        B(i, c) = b;
        
        w = F0;
        w = subs(w, unknown, unknown0);
        w = subs(w,L,L0(c)');
        w = eval(w);
        W(i) = w;
    end
        
        
   
    W =W+B*(L00-L0);
    dx = vpa(-pinv(A'*pinv(B*Qy*B')*A)*A'*pinv(B*Qy*B')*W);
    e0=e1;
    e1=vpa(Qy*B'*pinv(B*Qy*B')*(eye(n)-A*pinv(A'*pinv(B*Qy*B')*A)*A'*pinv(B*Qy*B'))*W);
    
    unknown0 = vpa(unknown0 + dx');
%     unknown0=[wrapTo2Pi(unknown0(1:3)) , unknown0(3:4)]
    L0=L00-e1;
    disp('Next repeatition')
       norm(dx)
end
unknown0=double(unknown0);
omega2=wrapTo2Pi(unknown0(1));phi2=wrapTo2Pi(unknown0(2));kapa2=wrapTo2Pi(unknown0(3));by=unknown0(4);bz=unknown0(5);
% disp('     omega2                phi2                    kapa2                  by                 bz')
unknown0=[unknown0, bx];


xxx=[XY_L(:,1)-xp XY_L(:,2)-yp -f*ones(n,1)];
x_1=xxx(:,1);
y_1=xxx(:,2);
z_1=xxx(:,3);
XXX=[XY_R(:,1)-xp XY_R(:,2)-yp -f*ones(n,1)]*eval(R');
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
parameters=['omega2 phi2 kapa2 by bz bx']
R=eval(R);
r1=((-f*ones(n,1).*(Xm./Zm)+xp*ones(n,1)-XY_L(:,1)).^2+(-f*ones(n,1).*(Ym./Zm)+yp*ones(n,1)-XY_L(:,2)).^2).^0.5;
r2=(((-f*ones(n,1).*(R(:,1)'*[[Xm-bx]';[Ym-by]';[Zm-bz]'])'./(R(:,3)'*[[Xm-bx]';[Ym-by]';[Zm-bz]'])')+xp*ones(n,1)-XY_R(:,1)).^2 + ((-f*ones(n,1).*(R(:,2)'*[[Xm-bx]';[Ym-by]';[Zm-bz]'])'./(R(:,3)'*[[Xm-bx]';[Ym-by]';[Zm-bz]'])')+yp*ones(n,1)-XY_R(:,2)).^2).^0.5;
pixel_error=(r1+r2)./2;
pixel_error=pixel_error./Pixel_Size;
end