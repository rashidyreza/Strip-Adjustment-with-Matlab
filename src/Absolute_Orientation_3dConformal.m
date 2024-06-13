function [EQ,RMSE,error]= Absolute_Orientation_3dConformal(X_Y_Z_1_Model,X_Y_Z_2_Object)
% X_Y_Z_1_Model=X_Y_Z_1_Model'
% X_Y_Z_2_Object=X_Y_Z_2_Object'
%=====================================
% X_Y_Z_2_Object=(Landa*Rk*Rph*Ro*X_Y_Z_1_Model)+T;
%=====================================
n=size(X_Y_Z_1_Model,2);
%=====================================
% syms h11 h12 h13 Tx h21 h22 h23 Ty h31 h32 h33 Tz
% H=[h11 h12 h13; h21 h22 h23; h31 h32 h33];
% UnKnown=[h11 h12 h13 Tx h21 h22 h23 Ty h31 h32 h33 Tz];
% q1=sum(H(:,1).*H(:,1))-sum(H(:,2).*H(:,2));
% q2=sum(H(:,1).*H(:,1))-sum(H(:,3).*H(:,3));
% q3=sum(H(:,1).*H(:,2));
% q4=sum(H(:,1).*H(:,3));
% q5=sum(H(:,2).*H(:,3));
% d0=[q1;q2;q3;q4;q5];
% d=jacobian(d0,UnKnown);
%=====================================
A=[X_Y_Z_1_Model',ones(n,1),zeros(n,8);
zeros(n,4),X_Y_Z_1_Model',ones(n,1),zeros(n,4);
zeros(n,8),X_Y_Z_1_Model',ones(n,1)];
L=[X_Y_Z_2_Object(1,:)';X_Y_Z_2_Object(2,:)';X_Y_Z_2_Object(3,:)'];
x=inv(A'*A)*A'*L;
%=====================================
% dx=1;
% while norm(dx)>1e-10
% x0=x;
% D=double(subs(d,UnKnown,x0'));
% x=inv(A'*A)*A'*(L)+inv(A'*A)*D'*inv(D*inv(A'*A)*D')*(-double(subs(d0,UnKnown,x0'))+D*x0-D*inv(A'*A)*A'*(L));
% dx=x-x0;
% end
abs_Matrix=reshape(x,4,3)';
% % % Landa0=sqrt(sum(abs_Matrix(1:3,1).^2));
% % % phi0=asin(abs_Matrix(3,1)/Landa0);
% % % omega0=asin(-abs_Matrix(3,2)*sec(phi0)/Landa0);
% % % kapa0=asin(-abs_Matrix(2,1)*sec(phi0)/Landa0);
% % % Tx0=abs_Matrix(1,4);
% % % Ty0=abs_Matrix(2,4);
% % % Tz0=abs_Matrix(3,4);
% % % syms kapa phi omega Landa Tx Ty Tz
% % % Rk=[cos(kapa) sin(kapa) 0; -sin(kapa) cos(kapa) 0; 0 0 1];
% % % Rph=[cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
% % % Romg=[1 0 0; 0 cos(omega) sin(omega); 0 -sin(omega) cos(omega)];
% % % H=Landa*Rk*Rph*Romg;
% % % EQ0=reshape(H*X_Y_Z_1_Model +[Tx;Ty;Tz]*ones(1,n),[],1)-X_Y_Z_2_Object(:);
% % % unknown = [kapa phi omega Landa Tx Ty Tz];
% % %  unknown0 = [kapa0 phi0 omega0 Landa0 Tx0 Ty0 Tz0];
% % % A0 = jacobian(EQ0,unknown);
% % % dx=1;
% % % while (norm(dx)> 1e-10)
% % %     A = eval(subs(A0,unknown,unknown0));
% % %     W = eval(subs(EQ0,unknown,unknown0));
% % %     dx = -inv(A'*A)*A'*W;
% % %     unknown0 = unknown0 + dx';
% % %     norm(dx)
% % % end
% unknown0-[kapa0 phi0 omega0 Landa0 Tx0 Ty0 Tz0]
% % kapa=unknown0(1);phi=unknown0(2);omega=unknown0(3);Landa=unknown0(4);Tx=unknown0(5);Ty=unknown0(6);Tz=unknown0(7);
% % % %-------------------------------------------------
% % % Rk=[cos(kapa) sin(kapa) 0; -sin(kapa) cos(kapa) 0; 0 0 1];
% % % Rph=[cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
% % % Romg=[1 0 0; 0 cos(omega) sin(omega); 0 -sin(omega) cos(omega)];
% % % EQ=[Landa*Rk*Rph*Romg,[Tx;Ty;Tz];0 0 0 1];
EQ=[abs_Matrix;0 0 0 1];
%RMSE
error=(EQ*[X_Y_Z_1_Model;ones(1,n)]-[X_Y_Z_2_Object;ones(1,n)]).^2;
error=sqrt(sum(error(1:3,:)))';
 RMSE=sqrt(sum(sum((EQ*[X_Y_Z_1_Model;ones(1,n)]-[X_Y_Z_2_Object;ones(1,n)]).^2)))/(3*n);
end