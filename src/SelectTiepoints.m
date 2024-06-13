function [Points_Im_1,Points_Im_2]= SelectTiepoints(I1,I2)

% I1=imread('B1_2.JPG');
% I2 = imread('B1_3.JPG');
II1=I1;
II2=I2;
I1 = rgb2gray(I1);
I2 = rgb2gray(I2);
points1 = detectHarrisFeatures(I1);
 points1=points1.selectStrongest(1500);
points2 = detectHarrisFeatures(I2);
  points2=points2.selectStrongest(1500);
[features1,valid_points1] = extractFeatures(I1,points1);
[features2,valid_points2] = extractFeatures(I2,points2);

indexPairs = matchFeatures(features1,features2);

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);


Points_Im_1=double(matchedPoints1.Location);
Points_Im_2=double(matchedPoints2.Location);
h=cpselect(II1,II2,Points_Im_1,Points_Im_2);
% [ZZZ,ZZZZ] = ...
%        cpselect(II1,II2,...
%                 Points_Im_1,Points_Im_2,...
%                 'Wait',true);
% menu('','Continue')


N=size(Points_Im_2,1);
X = java_array('java.lang.String', 3);
for i =1:N
X(i) = java.lang.String(char(string('P')+string(i)));
end
X(N+1)=java.lang.String(char(string('End')));
X(N+2)=java.lang.String(char(string('ALL')));
X(N+3)=java.lang.String(char(string('Remove some points')));
D = cell(X);
CHoice=0;
K=1;
while ( CHoice~=(N+1) && CHoice~=(N+2) ) && (CHoice~=(N+3))
CHoice=menu('Choose the points you want',D);
n(K)=CHoice;
K=K+1;
end
if CHoice==N+1
n(K-1)=[];
Points_Im_1=Points_Im_1(n,:);
Points_Im_2=Points_Im_2(n,:);
end

if CHoice==N+3
K=1;n=[];
while CHoice~=(N+1)
CHoice=menu('Remove the points you want',D(1:N+1));
n(K)=CHoice;
K=K+1;
end
n(K-1)=[];
Points_Im_1(n,:)=[];
Points_Im_2(n,:)=[];
end

close all
close(h)
showMatchedFeatures(II1,II2,[size(I1,2)/2 size(I1,1)/2],[size(I1,2)/2 size(I1,1)/2],'montage','Parent',axes);
CCC=menu('do you want to select more point','yes','no');
    while (CCC==1)
    'Celect1_in1'
    zoom
    pause
    centre1_in1=ginput(1);
    'Celect1_in2'
    zoom
    pause
    centre1_in2=ginput(1);
ccc=menu('do you want to add them','yes','no');
if ccc==1
Points_Im_1=[centre1_in1;Points_Im_1];
Points_Im_2=[[centre1_in2-[size(I2,2) 0]];Points_Im_2]; 
end
 CCC=menu('do you want to select more point','yes','no');
    end
[Points_Im_1,Points_Im_2] = ...
       cpselect(II1,II2,...
                Points_Im_1,Points_Im_2,...
                'Wait',true);
end