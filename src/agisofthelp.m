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
Points_Im_1_Pixel_Cs00=[x y]
Points_Im_2_Pixel_Cs00=[x1 y1]