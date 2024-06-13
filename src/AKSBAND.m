for i=7:1;%A.O for first image of each band
    COlumn_of_image=(Number_of_images_in_each_band*(j-1)+i)*2+4;% first image
    A=Tie_points_for_AO(:,COlumn_of_image:COlumn_of_image+1);
    aa=find(A(:,1)~=0);%rows id of gcp in first image
    B=Tie_points_for_AO(:,COlumn_of_image-2:COlumn_of_image-1);%second image
    bb=find(B(:,1)~=0);%rows id of gcp in second image
    C=Tie_points_for_AO(:,COlumn_of_image-4:COlumn_of_image-3);%third image
    cc=find(C(:,1)~=0);%rows id of gcp in third image
    dd=intersect(aa,bb);%rows id of gcp in first and second image
    dd=intersect(dd,cc);%rows id of gcp in first and second and third image
    XY_R=A(dd,:).*0.001;
    XY_L=B(dd,:).*0.001;
    unknown0=Band(j).Models(i+1).parameter;
%     XYZ_O_B=[unknown0(6) unknown0(4:5)]
       XYZ_O_B=[0 0 0];
   
    X_Y_Z_fixed_Model=Triangulation(unknown0,XY_L,XY_R,f,xp,yp,k1,k2,p1,p2); %Reference Model 1 i
% % %     X_Y_Z_fixed_Model=Absolute_Orientation_matrix(j).MAT(i+1).MATRIX*[X_Y_Z_fixed_Model';ones(1,size(X_Y_Z_fixed_Model,1))];
% % %     XYZ_O_B=Absolute_Orientation_matrix(j).MAT(i+1).MATRIX*[XYZ_O_B';1];
% % %     XYZ_O_B=XYZ_O_B(1:3)
% % %     X_Y_Z_fixed_Model=X_Y_Z_fixed_Model(1:3,:);
% % %     X_Y_Z_fixed_Model=X_Y_Z_fixed_Model';
    XY_R=B(dd,:).*0.001;
    XY_L=C(dd,:).*0.001;
    unknown0=Band(j).Models(i).parameter;
%      XYZ_O_M=[0 0 0];
     XYZ_O_M=[unknown0(6) unknown0(4:5)]
    X_Y_Z_Moving_Model=Triangulation(unknown0,XY_L,XY_R,f,xp,yp,k1,k2,p1,p2); %Moving Model 2 i+1
%    [EQ,RMSE,error]=Absolute_Orientation_3dConformal(X_Y_Z_Moving_Model',X_Y_Z_fixed_Model')
   
   [EQ,RMSE,error]=Absolute_Orientation_3dConformal_M43(X_Y_Z_Moving_Model',XYZ_O_M',X_Y_Z_fixed_Model',XYZ_O_B')

   
   
    sk=0.6
        while (max(error)> allow) && (size(find(error<sk*max(error)),1)>10)
        X_Y_Z_fixed_Model=X_Y_Z_fixed_Model(find(error<sk*max(error)),:);
        X_Y_Z_Moving_Model=X_Y_Z_Moving_Model(find(error<sk*max(error)),:);
                  [EQ,RMSE,error]=Absolute_Orientation_3dConformal_M43(X_Y_Z_Moving_Model',XYZ_O_M',X_Y_Z_fixed_Model',XYZ_O_B')
          [EQ,RMSE,error]=Absolute_Orientation_3dConformal(X_Y_Z_Moving_Model',X_Y_Z_fixed_Model')
   end
field1='MATRIX';
field2='RMSE';
field3='error';
% saving absolute orientation result
Result=struct(field1,EQ,field2,RMSE,field3,error);
%---------Number of band----------------
 Absolute_Orientation_matrix(j).MAT(i)=Result;
end
for i=7:1
 Absolute_Orientation_matrix(j).MAT(i).MATRIX=Absolute_Orientation_matrix(j).MAT(i+1).MATRIX* Absolute_Orientation_matrix(j).MAT(i).MATRIX
end