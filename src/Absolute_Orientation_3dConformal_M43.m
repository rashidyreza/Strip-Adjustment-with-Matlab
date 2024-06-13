function [EQ,RMSE,error]= Absolute_Orientation_3dConformal_M43(X_Y_Z_1_Model,XYZ_O_M,X_Y_Z_2_Object,XYZ_O_B)

% XYZ_O_M=XYZ_O_M'
% X_Y_Z_1_Model=X_Y_Z_Moving_Model'
% X_Y_Z_2_Object=X_Y_Z_fixed_Model'
%=====================================
% X_Y_Z_2_Object=(Landa*Rk*Rph*Ro*X_Y_Z_1_Model)+T;
%=====================================
X_Y_Z_1_Model_0=X_Y_Z_1_Model;
XYZ_O_M_0=XYZ_O_M;
EQ=eye(4);
n=size(X_Y_Z_1_Model,2);
%=====================================
MAT=5*eye(4);
for iii=1:6
% while norm(MAT-eye(4))>1e-2
% h=1;
% while h~=4
%     h=h+1;
G_Model=X_Y_Z_1_Model(1:2,:)-mean(X_Y_Z_1_Model(1:2,:)')'*ones(1,n);
G_Object=X_Y_Z_2_Object(1:2,:)-mean(X_Y_Z_2_Object(1:2,:)')'*ones(1,n);
% % G_Model=X_Y_Z_1_Model(1:2,:);
% % G_Object=X_Y_Z_2_Object(1:2,:);
A=[];L=[];
A(1:2:2*n,:)=[ G_Model(1,:)' G_Model(2,:)' ,ones(n,1) ,zeros(n,1)];
A(2:2:2*n,:)=[G_Model(2,:)' -G_Model(1,:)' ,zeros(n,1) ,ones(n,1)];
L(1:2:2*n,:)=G_Object(1,:)';
L(2:2:2*n,:)=G_Object(2,:)';
N=A'*A;u=A'*L;
a=u(1)/N(1,1);
b=u(2)/N(2,2);
landa=sqrt(a^2+b^2);
kapa=atan2(b,a);
Cx= mean(X_Y_Z_2_Object(1,:)-a.*X_Y_Z_1_Model(1,:)-b.*X_Y_Z_1_Model(2,:));
Cy= mean(X_Y_Z_2_Object(2,:)-a.*X_Y_Z_1_Model(2,:)+b.*X_Y_Z_1_Model(1,:));
% % unknown_M4=inv(N)*u;
% % a=unknown_M4(1);b=unknown_M4(2);Cx=unknown_M4(3);Cy=unknown_M4(4); landa=sqrt(a^2+b^2);
% %  kapa=atan(b/a);
%=====================================
X_Y_Z_second_Model=[a b 0;-b a 0;0 0 landa]*X_Y_Z_1_Model+[Cx;Cy;0]*ones(1,n);
XYZ_O_M=[a b 0;-b a 0;0 0 landa]*XYZ_O_M+[Cx;Cy;0];
%=====================================W_Phi_Cz
A=[];
L=[];
A=[X_Y_Z_second_Model(2,:)' -X_Y_Z_second_Model(1,:)' ,ones(n,1)];
L=[X_Y_Z_2_Object(3,:)-X_Y_Z_second_Model(3,:)]';
% L=[L;XYZ_O_B-XYZ_O_M];
% A=[A;0 XYZ_O_M(3) 0;-XYZ_O_M(3) 0 0;XYZ_O_M(2) -XYZ_O_M(1) 1];
Unknown_M3=inv(A'*A)*A'*L;
W=Unknown_M3(1);
phi=Unknown_M3(2);
Cz=Unknown_M3(3);


R_kapa = [cos(kapa) sin(kapa) 0 ; -sin(kapa) cos(kapa) 0 ; 0 0 1];
R_phi = [cos(phi) 0 sin(phi) ; 0 1 0 ; -sin(phi) 0 cos(phi)];
R_omega = [1 0 0 ; 0 cos(W) -sin(W) ; 0 sin(W) cos(W)];
R = R_phi*R_omega;
%% 
% syms Cz phi W
% R_phi = [cos(phi) 0 sin(phi) ; 0 1 0 ; -sin(phi) 0 cos(phi)];
% R_omega = [1 0 0 ; 0 cos(W) -sin(W) ; 0 sin(W) cos(W)];
% R = R_phi*R_omega;
% F=X_Y_Z_2_Object(3,:)-R(3,:)*X_Y_Z_second_Model-Cz*ones(1,n);
% F=reshape(F,[],1);
% unknown = [W phi Cz];
% unknown0 = [0 0 1];
% A0 = jacobian(F,unknown);
% dx=1;
% while (norm(dx)> 1e-10)
%     A = eval(subs(A0,unknown,unknown0));
%     df = eval(subs(F,unknown,unknown0));
%     dx =-inv(A'*A)*A'*df;
%     unknown0 = unknown0 + dx';
%     norm(dx)
% end
% W=unknown0(1);phi=unknown0(2);Cz=unknown0(3);
% R_phi = [cos(phi) 0 sin(phi) ; 0 1 0 ; -sin(phi) 0 cos(phi)];
% R_omega = [1 0 0 ; 0 cos(W) -sin(W) ; 0 sin(W) cos(W)];
% R = R_phi*R_omega;
% R_kapa = [cos(kapa) sin(kapa) 0 ; -sin(kapa) cos(kapa) 0 ; 0 0 1];
%==========================Apply M3==========
% X_Y_Z_1_Model=[1 0 phi;0 1 -W;-phi W 1]*X_Y_Z_second_Model+[0;0;Cz]*ones(1,n);
% XYZ_O_M=[1 0 phi;0 1 -W;-phi W 1]*XYZ_O_M+[0;0;Cz];
% MAT=[1 0 phi 0;0 1 -W 0;-phi W 1 Cz;0 0 0 1]*[a b 0 Cx;-b a 0 Cy;0 0 landa 0;0 0 0 1]
% EQ=MAT*EQ;
X_Y_Z_1_Model=R*X_Y_Z_second_Model+[0;0;Cz]*ones(1,n);
XYZ_O_M=R*XYZ_O_M+[0;0;Cz];
MAT=[[R;0 0 0],[0 ;0 ;Cz ;1]]*[[landa*R_kapa;0 0 0],[Cx ;Cy ;0 ;1]];
EQ=MAT*EQ;
end
phi
Cz
W
kapa
Cx
Cy
landa
%RMSE
error=(EQ*[X_Y_Z_1_Model_0;ones(1,n)]-[X_Y_Z_2_Object;ones(1,n)]).^2;
error=sqrt(sum(error(1:3,:)))';
 RMSE=sqrt(sum(sum((EQ*[X_Y_Z_1_Model_0;ones(1,n)]-[X_Y_Z_2_Object;ones(1,n)]).^2)))/(3*n);
%  error=X_Y_Z_1_Model-X_Y_Z_2_Object
end