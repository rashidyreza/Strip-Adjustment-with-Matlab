[n1,n2]=uigetfile('*.txt','Select Export matches');
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
XYZ_GCP=dlmread('GCP_UTM.txt')
XYZ_GCP(:,2:4)=XYZ_GCP(:,2:4)- [606850.*ones(size(XYZ_GCP,1),1)  3486400.*ones(size(XYZ_GCP,1),1) 1600.*ones(size(XYZ_GCP,1),1)]
for i=1:size(Tie_points_for_AO,1)
    Tie_points_for_AO(i,Number_of_Bands*Number_of_images_in_each_band*2+2:Number_of_Bands*Number_of_images_in_each_band*2+4)=XYZ_GCP(find(XYZ_GCP(:,1)==Tie_points_for_AO(i,1)),2:4);
    Tie_points_for_AO_pixel(i,Number_of_Bands*Number_of_images_in_each_band*2+2:Number_of_Bands*Number_of_images_in_each_band*2+4)=XYZ_GCP(find(XYZ_GCP(:,1)==Tie_points_for_AO(i,1)),2:4);
end



