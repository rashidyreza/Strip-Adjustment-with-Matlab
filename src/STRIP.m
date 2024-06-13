clc;clear;close all;
Pixel_Size=0.00156192*0.001;
f=2312.78104*Pixel_Size;
xp=0;
yp=0;
k1=8.756635940e-004;
k2=-8.843480548e-005;
p1=0;
p2=0;
rows=3000;
colums=4000;
Number_of_Bands=3;
Number_of_images_in_each_band=10;
KK=1;
Models=[];
for j=1:Number_of_Bands
for i=1:Number_of_images_in_each_band-1
Im1=string('B')+string(j)+string('_')+string(i)+string('.JPG');
Im2=string('B')+string(j)+string('_')+string(i+1)+string('.JPG');
Im1=char(Im1);
Im2=char(Im2);
I1=imread(Im1);
I2=imread(Im2);
CHoice=menu('Select matches','Matlab_matches','get help from agisoft');
switch CHoice
    case 1
[Points_Im_1_Pixel_Cs00,Points_Im_2_Pixel_Cs00]= SelectTiepoints(I1,I2);
    case 2
[n1,n2]=uigetfile('*.txt','Select Export matches');
Export_matches=strcat(n2,n1);
Export_matches=readtable(Export_matches);
A=Export_matches(find(Export_matches.Var1==string(Im1)),2:4);
B=Export_matches(find(Export_matches.Var1==string(Im2)),2:4);
C=intersect(A(:,1),B(:,1));
rr=A(find(ismember(A(:,1),C)==1),:);
rr1=B(find(ismember(B(:,1),C)==1),:);
rr=table2array(rr(:,2:3));
rr=rr/(Pixel_Size*1000);
x=rr(:,1)+2000+0.5;
y=-rr(:,2)+1500+0.5;
rr1=table2array(rr1(:,2:3));
rr1=rr1./(Pixel_Size*1000);
x1=rr1(:,1)+colums/2+0.5;
y1=-rr1(:,2)+rows/2+0.5;
choice=menu('Do you want to see matches','yes','no')
switch choice
    case 1
imshow(I1)
hold on
plot(x,y,'*')
%---------------------
figure;
imshow(I2)
hold on
plot(x1,y1,'*')
end
Points_Im_1_Pixel_Cs00=[x y];
Points_Im_2_Pixel_Cs00=[x1 y1];
end
Points_Im_1_Pixel_Cs=Points_Im_1_Pixel_Cs00;
Points_Im_2_Pixel_Cs=Points_Im_2_Pixel_Cs00;
%----------------Image Coordinate system-----------------------------
Points_Im_1_Image_Cs= [ Points_Im_1_Pixel_Cs00(:,1)-colums/2  rows/2-Points_Im_1_Pixel_Cs00(:,2) ].*Pixel_Size;
Points_Im_2_Image_Cs= [ Points_Im_2_Pixel_Cs00(:,1)-colums/2  rows/2-Points_Im_2_Pixel_Cs00(:,2) ].*Pixel_Size;
%----------------Relative_Orientation---------Model Coordinate system-------------------------

[X_Y_Z_Model,var_py,py,pixel_error,parameters,unknown0]=Relative_Orientation_Triangulation2(Points_Im_1_Image_Cs,Points_Im_2_Image_Cs,f,xp,yp,Pixel_Size,k1,k2,p1,p2)
disp(std(pixel_error))


%--------------------------------Fix-----pixel_error--------------------------------------
allow_Pe=0.4;
% % % error_Points_Im_1_Pixel_Cs=Points_Im_1_Pixel_Cs(find(pixel_error>std(pixel_error)),:);
% % % error_Points_Im_2_Pixel_Cs=Points_Im_2_Pixel_Cs(find(pixel_error>std(pixel_error)),:);
% % % if size(error_Points_Im_1_Pixel_Cs,1)==0
% % %     error_Points_Im_1_Pixel_Cs=[0 0];
% % %     error_Points_Im_2_Pixel_Cs=[0 0];
% % % end
% % % Points_Im_1_Pixel_Cs(find(pixel_error>std(pixel_error)),:)=[];
% % % Points_Im_2_Pixel_Cs(find(pixel_error>std(pixel_error)),:)=[];
% % % [Fixed_Points_Im_1_Pixel_Cs,Fixed_Points_Im_2_Pixel_Cs] = ...
% % %        cpselect(I1,I2,...
% % %                 error_Points_Im_1_Pixel_Cs,error_Points_Im_2_Pixel_Cs,...
% % %                 'Wait',true);
% % % Points_Im_1_Pixel_Cs=[Points_Im_1_Pixel_Cs;Fixed_Points_Im_1_Pixel_Cs];
% % % Points_Im_2_Pixel_Cs=[Points_Im_2_Pixel_Cs;Fixed_Points_Im_2_Pixel_Cs];
% % % Points_Im_1_Image_Cs= [ Points_Im_1_Pixel_Cs(:,1)-colums/2  rows/2-Points_Im_1_Pixel_Cs(:,2) ].*Pixel_Size;
% % % Points_Im_2_Image_Cs= [ Points_Im_2_Pixel_Cs(:,1)-colums/2  rows/2-Points_Im_2_Pixel_Cs(:,2) ].*Pixel_Size;
% % % [X_Y_Z_Model,var_py,py,pixel_error,parameters,unknown0]=Relative_Orientation_Triangulation2(Points_Im_1_Image_Cs,Points_Im_2_Image_Cs,f,xp,yp,Pixel_Size,k1,k2,p1,p2)
% % % disp(std(pixel_error))
%-------------------------------------------------
while std(pixel_error)> 0.5*allow_Pe && (size(Points_Im_1_Image_Cs(find(pixel_error<std(pixel_error))),1)>=15)
while (max(pixel_error) > allow_Pe)&& (size(find(pixel_error)<std(pixel_error),1)>=15)
Points_Im_1_Image_Cs=Points_Im_1_Image_Cs(find(pixel_error<2.5*std(pixel_error)),:);
Points_Im_2_Image_Cs=Points_Im_2_Image_Cs(find(pixel_error<2.5*std(pixel_error)),:);
Points_Im_1_Pixel_Cs=Points_Im_1_Pixel_Cs(find(pixel_error<2.5*std(pixel_error)),:);
Points_Im_2_Pixel_Cs=Points_Im_2_Pixel_Cs(find(pixel_error<2.5*std(pixel_error)),:);
[X_Y_Z_Model,var_py,py,pixel_error,parameters,unknown0]=Relative_Orientation_Triangulation2(Points_Im_1_Image_Cs,Points_Im_2_Image_Cs,f,xp,yp,Pixel_Size,k1,k2,p1,p2)
disp(std(pixel_error))
end
end
if (size(Points_Im_1_Image_Cs(find(pixel_error>0.5*allow_Pe),:),1)>0) && (size(Points_Im_1_Image_Cs(find(pixel_error<0.5*allow_Pe)),1)>10)
Points_Im_1_Image_Cs=Points_Im_1_Image_Cs(find(pixel_error<0.5*allow_Pe),:);
Points_Im_2_Image_Cs=Points_Im_2_Image_Cs(find(pixel_error<0.5*allow_Pe),:);
Points_Im_1_Pixel_Cs=Points_Im_1_Pixel_Cs(find(pixel_error<0.5*allow_Pe),:);
Points_Im_2_Pixel_Cs=Points_Im_2_Pixel_Cs(find(pixel_error<0.5*allow_Pe),:);
[X_Y_Z_Model,var_py,py,pixel_error,parameters,unknown0]=Relative_Orientation_Triangulation2(Points_Im_1_Image_Cs,Points_Im_2_Image_Cs,f,xp,yp,Pixel_Size,k1,k2,p1,p2)
end
%---------------------------------
% % % h=cpselect(I1,I2,Points_Im_1_Pixel_Cs,Points_Im_2_Pixel_Cs);
% % % close(h);
%-------------------------------

parameter = unknown0;
ModelCs=X_Y_Z_Model;
Im1_ps_Cs=Points_Im_1_Pixel_Cs;
Im2_ps_Cs=Points_Im_2_Pixel_Cs;
Im1_Im_Cs=Points_Im_1_Image_Cs;
Im2_Im_Cs=Points_Im_2_Image_Cs;
%-------------------------------

field1='parameter';
field2='ModelCs';
field3='Im1_ps_Cs';
field4='Im2_ps_Cs';
field5='Im1_Im_Cs';
field6='Im2_Im_Cs';
field7='Points_Im_1_Pixel_Cs00';
field8='Points_Im_2_Pixel_Cs00';
field9='pixel_error';
% saving relative orientation parameter
MODELpa=struct(field1,unknown0,field2,X_Y_Z_Model,field3,Points_Im_1_Pixel_Cs,field4,Points_Im_2_Pixel_Cs,field5,Points_Im_1_Image_Cs,field6,Points_Im_2_Image_Cs,field7,Points_Im_1_Pixel_Cs00,field8,Points_Im_2_Pixel_Cs00,field9,pixel_error);
%---------Number of band----------------
 Band(j).Models(i)=MODELpa;
end
end
[n1,n2]=uigetfile('*.txt','GCP_ImageCs.txt');
T=strcat(n2,n1);
T=readtable(T);
GCPPoint_identification=zeros(size(unique(T.Num),1),Number_of_Bands*Number_of_images_in_each_band*2+1);
GCPPoint_identification(:,1)=unique(T.Num);
GCPPoint_identification_pixel=GCPPoint_identification;
ALL_point_IDs=T.Num;
for j=1:Number_of_Bands
for i=1:Number_of_images_in_each_band
COlumn_of_image=(Number_of_images_in_each_band*(j-1)+i)*2;
Im=string('B')+string(j)+string('_')+string(i)+string('.JPG');
Point_Location=find(T.image==Im);%Point_Location
Points_ID=ALL_point_IDs(find(T.image==Im));%Points_ID
for kk=1:size(Points_ID,1);
   satr=find(GCPPoint_identification(:,1)==Points_ID(kk));
   GCPPoint_identification(satr,COlumn_of_image)=T.x(Point_Location(kk))+Pixel_Size/2;
   GCPPoint_identification(satr,COlumn_of_image+1)=T.y(Point_Location(kk))-Pixel_Size/2;
   GCPPoint_identification_pixel(satr,COlumn_of_image)=T.x(Point_Location(kk))/(Pixel_Size*1000)+2000.5;
   GCPPoint_identification_pixel(satr,COlumn_of_image+1)=T.y(Point_Location(kk))/(Pixel_Size*1000)+1500.5;
end
end
end
XYZ_GCP=dlmread('GCP_UTM.txt');
XYZ_GCP(:,2:4)=XYZ_GCP(:,2:4)- [606850.*ones(size(XYZ_GCP,1),1)  3486400.*ones(size(XYZ_GCP,1),1) 1600.*ones(size(XYZ_GCP,1),1)];
for i=1:size(GCPPoint_identification,1);
    GCPPoint_identification(i,Number_of_Bands*Number_of_images_in_each_band*2+2:Number_of_Bands*Number_of_images_in_each_band*2+4)=XYZ_GCP(find(XYZ_GCP(:,1)==GCPPoint_identification(i,1)),2:4);
    GCPPoint_identification_pixel(i,Number_of_Bands*Number_of_images_in_each_band*2+2:Number_of_Bands*Number_of_images_in_each_band*2+4)=XYZ_GCP(find(XYZ_GCP(:,1)==GCPPoint_identification(i,1)),2:4);
end
%%--------------------Tie points for Absolute Orientation------------------
[n1,n2]=uigetfile('*.txt','Select Export matches(Tie_points_for_AO)');
T=strcat(n2,n1);
T=readtable(T);
Tie_points_for_AO=zeros(size(unique(T.Num),1),Number_of_Bands*Number_of_images_in_each_band*2+1);
Tie_points_for_AO(:,1)=sort(unique(T.Num));
Tie_points_for_AO_pixel=Tie_points_for_AO;
ALL_point_IDs=T.Num;
for j=1:Number_of_Bands
for i=1:Number_of_images_in_each_band
COlumn_of_image=(Number_of_images_in_each_band*(j-1)+i)*2;
Im=string('B')+string(j)+string('_')+string(i)+string('.JPG');
Point_Location=find(T.image==Im);%Point_Location
Points_ID=ALL_point_IDs(find(T.image==Im));%Points_ID
for kk=1:size(Points_ID,1)
   satr=find(Tie_points_for_AO(:,1)==Points_ID(kk));
   Tie_points_for_AO(satr,COlumn_of_image)=T.x(Point_Location(kk))+Pixel_Size/2;
   Tie_points_for_AO(satr,COlumn_of_image+1)=T.y(Point_Location(kk))-Pixel_Size/2;
   Tie_points_for_AO_pixel(satr,COlumn_of_image)=T.x(Point_Location(kk))/(Pixel_Size*1000)+2000.5;
   Tie_points_for_AO_pixel(satr,COlumn_of_image+1)=T.y(Point_Location(kk))/(Pixel_Size*1000)+1500.5;
end
end
end


GCP=GCPPoint_identification(:,62:63);
GCP=GCP';
teta=2*pi+(atan2(607194.92-607378.04,3487026.184-3487235.911))-pi/2;
teta=-teta;
GCP=[cos(teta) sin(teta);-sin(teta) cos(teta)]*GCP;
GCP=GCP';
GCPPoint_identification(:,62:63)=GCP;


% % point_1_89_Ro=[cos(teta) sin(teta);-sin(teta) cos(teta)]*[607194.92 607378.04;3487026.184 3487235.911];
% % atan2(point_1_89_Ro(1,1)-point_1_89_Ro(1,2),point_1_89_Ro(2,1)-point_1_89_Ro(2,2))
% % [cos(teta) sin(-teta);sin(teta) cos(teta)]*[-3026038.10145938,-3026316.52263122;-1836063.9115487,-1836063.9115487]

% ---------------Absolute_Orientation_3dConformal----------------------------------
i=1;%A.O for first image of each band
allow=0.2;
for j=1:Number_of_Bands % first model in each band
    COlumn_of_image=(Number_of_images_in_each_band*(j-1)+i)*2;% first image
    A=GCPPoint_identification(:,COlumn_of_image:COlumn_of_image+1); % x y imL
    aa=find(A(:,1)~=0);%rows id of gcp in left image
    B=GCPPoint_identification(:,COlumn_of_image+2:COlumn_of_image+3); %x y imR
    bb=find(B(:,1)~=0);%rows id of gcp in right image
    cc=intersect(aa,bb);%rows id of gcp in both image
    X_Y_Z_2_Object=GCPPoint_identification(cc,Number_of_Bands*Number_of_images_in_each_band*2+2:Number_of_Bands*Number_of_images_in_each_band*2+4);
    XY_L=A(cc,:).*0.001;
    XY_R=B(cc,:).*0.001;
    unknown0=Band(j).Models(i).parameter;
    [X_Y_Z_1_Model]=Triangulation(unknown0,XY_L,XY_R,f,xp,yp,k1,k2,p1,p2);
    [EQ,RMSE,error]=Absolute_Orientation_3dConformal(X_Y_Z_1_Model',X_Y_Z_2_Object')
    while (max(error)> allow) && (size(find(error<max(error)),1)>6)
        X_Y_Z_1_Model=X_Y_Z_1_Model(find(error<max(error)),:);
        X_Y_Z_2_Object=X_Y_Z_2_Object(find(error<max(error)),:);
         [EQ,RMSE,error]=Absolute_Orientation_3dConformal(X_Y_Z_1_Model',X_Y_Z_2_Object')
    end
field1='MATRIX';
field2='RMSE';
field3='error';
% saving absolute orientation result
Result=struct(field1,EQ,field2,RMSE,field3,error);
%---------Number of band----------------
 Absolute_Orientation_matrix(j).MAT(i)=Result;
end
% ---------------Absolute_Orientation_3dConformal----------------------------------
allow=0.2;
for j=1:Number_of_Bands % first model in each band
for i=1:Number_of_images_in_each_band-2;%A.O for first image of each band
    COlumn_of_image=(Number_of_images_in_each_band*(j-1)+i)*2;% first image
    A=Tie_points_for_AO(:,COlumn_of_image:COlumn_of_image+1);
    aa=find(A(:,1)~=0);%rows id of gcp in first image
    B=Tie_points_for_AO(:,COlumn_of_image+2:COlumn_of_image+3);%second image
    bb=find(B(:,1)~=0);%rows id of gcp in second image
    C=Tie_points_for_AO(:,COlumn_of_image+4:COlumn_of_image+5);%third image
    cc=find(C(:,1)~=0);%rows id of gcp in third image
    dd=intersect(aa,bb);%rows id of gcp in first and second image
    dd=intersect(dd,cc);%rows id of gcp in first and second and third image
    XY_L=A(dd,:).*0.001;
    XY_R=B(dd,:).*0.001;
    unknown0=Band(j).Models(i).parameter;
    X_Y_Z_fixed_Model=Triangulation(unknown0,XY_L,XY_R,f,xp,yp,k1,k2,p1,p2); %Reference Model 1 i
    X_Y_Z_fixed_Model=Absolute_Orientation_matrix(j).MAT(i).MATRIX*[X_Y_Z_fixed_Model';ones(1,size(X_Y_Z_fixed_Model,1))];
    X_Y_Z_fixed_Model=X_Y_Z_fixed_Model(1:3,:);
    X_Y_Z_fixed_Model=X_Y_Z_fixed_Model';
    XY_L=B(dd,:).*0.001;
    XY_R=C(dd,:).*0.001;
    unknown0=Band(j).Models(i+1).parameter;
    X_Y_Z_Moving_Model=Triangulation(unknown0,XY_L,XY_R,f,xp,yp,k1,k2,p1,p2); %Moving Model 2 i+1
    [EQ,RMSE,error]=Absolute_Orientation_3dConformal(X_Y_Z_Moving_Model',X_Y_Z_fixed_Model')
        while (max(error)> allow) && (size(find(error<max(error)),1)>7)
        X_Y_Z_fixed_Model=X_Y_Z_fixed_Model(find(error<max(error)),:);
        X_Y_Z_Moving_Model=X_Y_Z_Moving_Model(find(error<max(error)),:);
         [EQ,RMSE,error]=Absolute_Orientation_3dConformal(X_Y_Z_Moving_Model',X_Y_Z_fixed_Model')
    end
field1='MATRIX';
field2='RMSE';
field3='error';
% saving absolute orientation result
Result=struct(field1,EQ,field2,RMSE,field3,error);
%---------Number of band----------------
 Absolute_Orientation_matrix(j).MAT(i+1)=Result;
end
end



% % % % load('Band_Are_ready_Absolute_Relative.mat')



Band_CS_identification=zeros(size(GCPPoint_identification,1),Number_of_Bands*3+1);
Band_CS_identification(:,1)=GCPPoint_identification(:,1);
for k=1:size(GCPPoint_identification,1)
    for j=1:Number_of_Bands
       A=GCPPoint_identification(k,Number_of_images_in_each_band*(j-1)*2+2:Number_of_images_in_each_band*(j)*2+1);
       AA= find(A(1:2:end)~=0);
        if (isempty(AA)~=1) && (size(AA,2)>1)
         AA=AA(find(([AA,0]-[0 AA+1])==0));
         if (isempty(AA)~=1)
         AA=AA(1); %Right
    XY_L=A(2*AA-3:2*AA-2).*0.001;
    XY_R=A(2*AA-1:2*AA).*0.001;
    unknown0=Band(j).Models(AA-1).parameter;
    X_Y_Z_fixed_Model=Triangulation(unknown0,XY_L,XY_R,f,xp,yp,k1,k2,p1,p2); %Reference Model 1 i
    X_Y_Z_fixed_Model=Absolute_Orientation_matrix(j).MAT(AA-1).MATRIX*[X_Y_Z_fixed_Model';ones(1,size(X_Y_Z_fixed_Model,1))];
    X_Y_Z_fixed_Model=X_Y_Z_fixed_Model(1:3,:);
    X_Y_Z_fixed_Model=X_Y_Z_fixed_Model';
    Band_CS_identification(k,3*j-1:3*j+1)=X_Y_Z_fixed_Model;
         end   
        end
    end
end
%%=======Strip Adj========
A=[];
B=[];
C=[];
I=-eye(3);
% gcp_1=[31;32;33;37;38;39;40];
gcp_1=[21;22;29;31;32;33;34;35;37;38;39;40;41;42;107]
% gcp_2=[31;32;33;37;38;39;40];
gcp_2=[21;22;28;29;30;31;32;33;34;35;37;38;39;40;41;107;108];
% gcp_3=[31;32;33;35;43;75;78];
gcp_3=[21;22;31;32;33;35;43;75;78];
tie_1_2=[17;18;24;27];
tie_1_3=[26;36;45;77];
%Strip Band 1
%gcp
for i=1:size(gcp_1,1)
    xyz=Band_CS_identification(find(Band_CS_identification(:,1)==gcp_1(i)),2:4);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    A(3*i-2:3*i,1:11)=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y];
    xyz=GCPPoint_identification(find(GCPPoint_identification(:,1)==gcp_1(i)),62:64);
    C(3*i-2:3*i,1)=xyz';
end
%tie 1 & 2
for i=1:size(tie_1_2,1)
    xyz=Band_CS_identification(find(Band_CS_identification(:,1)==tie_1_2(i)),2:4);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    A(3*size(gcp_1,1)+3*i-2:3*size(gcp_1,1)+3*i,1:11)=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y];
    B(3*size(gcp_1,1)+3*i-2:3*size(gcp_1,1)+3*i,3*i-2:3*i)=I;
end
% tie 1 & 3
for i=1:size(tie_1_3,1)
    xyz=Band_CS_identification(find(Band_CS_identification(:,1)==tie_1_3(i)),2:4);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    A(3*size(tie_1_2,1)+3*size(gcp_1,1)+3*i-2:3*size(tie_1_2,1)+3*size(gcp_1,1)+3*i,1:11)=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y];
    B(3*size(tie_1_2,1)+3*size(gcp_1,1)+3*i-2:3*size(tie_1_2,1)+3*size(gcp_1,1)+3*i,3*size(tie_1_2,1)+3*i-2:3*size(tie_1_2,1)+3*i)=I;
end
k=3*size(tie_1_3,1)+3*size(tie_1_2,1)+3*size(gcp_1,1);
%Strip Band 2
%gcp
for i=1:size(gcp_2,1)
    xyz=Band_CS_identification(find(Band_CS_identification(:,1)==gcp_2(i)),5:7);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    A(k+3*i-2:k+3*i,12:22)=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y];
    xyz=GCPPoint_identification(find(GCPPoint_identification(:,1)==gcp_2(i)),62:64);
    C(k+3*i-2:k+3*i,1)=xyz';
end
k=k+3*size(gcp_2,1);
%tie 1 & 2
for i=1:size(tie_1_2,1)
    xyz=Band_CS_identification(find(Band_CS_identification(:,1)==tie_1_2(i)),5:7);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    A(k+3*i-2:k+3*i,12:22)=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y];
    B(k+3*i-2:k+3*i,3*i-2:3*i)=I;
end
k=k+3*size(tie_1_2,1);
%Strip Band 3
%gcp
for i=1:size(gcp_3,1)
    xyz=Band_CS_identification(find(Band_CS_identification(:,1)==gcp_3(i)),8:10);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    A(k+3*i-2:k+3*i,23:33)=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y];
    xyz=GCPPoint_identification(find(GCPPoint_identification(:,1)==gcp_3(i)),62:64);
    C(k+3*i-2:k+3*i,1)=xyz';
end
k=k+3*size(gcp_3,1);
%tie 1 & 3
for i=1:size(tie_1_3,1)
    xyz=Band_CS_identification(find(Band_CS_identification(:,1)==tie_1_3(i)),8:10);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    A(k+3*i-2:k+3*i,23:33)=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y];
    B(k+3*i-2:k+3*i,3*size(tie_1_2,1)+3*i-2:3*size(tie_1_2,1)+3*i)=I;
    C(k+3*i-2:k+3*i,1)=[0,0,0]';
end
%@@@@@@@@@@@@@@@@ESTIMATION@@@@@@@@@@@@@@
N11=A'*A;
N12=A'*B;
F1=A'*C;
N21=B'*A;
N22=B'*B;
unknown=inv(N11-N21'*inv(N22)*N21)*(F1);
a=[A B]\C;
%--------------------CHECKING--------------------
point=371;
choice=menu('select your band','band 1','band 2','band 3');
switch choice
    case 1
b=1:11;c=2:4
    case 2
b=12:22;c=5:7
    case 3
b=23:33;c=8:10
end
  xyz=Band_CS_identification(find(Band_CS_identification(:,1)==point),c);
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
CHECK=[400,x,x^2,0,-y,-2*x*y,0,z,2*x*z,0,0;0,y,2*x*y,400,x,x^2,0,0,0,-z,-2*x*z;0,z,2*x*z,0,0,0,400,x,x^2,y,2*x*y]*unknown(b)
 xyz=GCPPoint_identification(find(GCPPoint_identification(:,1)==point),62:64);
CHECK-xyz'
%-------------------------------------------------------



ANS=Band_CS_identification(:,2:end)-[GCPPoint_identification(:,62:64) GCPPoint_identification(:,62:64) GCPPoint_identification(:,62:64)];
ANS=ANS.*(Band_CS_identification(:,2:end)~=0);
% teta=abs(atan2(607194.92-607378.04,3487026.184-3487235.911));




