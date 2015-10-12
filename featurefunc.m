
function [solid,axisarea,Con,Ene,Homo,Cor,convexx]  = featurefunc(S,maxArea1,BW2,im2)


  
idx1 = find([S.Area] >= maxArea1);
hull=S(idx1).ConvexHull;
solid=S(idx1).Solidity;
convexx=S(idx1).ConvexImage;
extrem=S(idx1).Extrema;%[top-left top-right right-top right-bottom bottom-right bottom-left left-bottom left-top]


x1=extrem(1,1);y1=extrem(1,2);
x2=extrem(2,1);y2=extrem(2,2);
x3=extrem(3,1);y3=extrem(3,2);
x4=extrem(4,1);y4=extrem(4,2);
x5=extrem(5,1);y5=extrem(5,2);
x6=extrem(6,1);y6=extrem(6,2);
x7=extrem(7,1);y7=extrem(7,2);
x8=extrem(8,1);y8=extrem(8,2);


if(x7<x8)
    x=x7;
else
    x=x8;
end

if(y1<y2)
    y=y1;
else
    y=y2;
end

if(x3>x4)
    xa=x3;
else
    xa=x4;
end

if(y5>y6)
    ya=y5;
else
    ya=y6;
end


rect=[x-10,y-10,xa-x+10,ya-y+10];
I=imcrop(BW2,rect);
% figure,imshow((I));
%imlabel('croped');m 
% vector[xmin ymin width height]
Imag=imresize(I,[50,50]);
% figure,imshow((Imag));
Imag1=imfill(Imag,'holes');
%figure,imshow(uint8(Imag1));

convexx=imresize(convexx,[200,100]);
con = regionprops(convexx,'all');
maxArea12 = max([con.Area]);
maxArea12=maxArea12-1;
idx1 = find([con.Area] >= maxArea12);
major=con(idx1).MajorAxisLength;


% figure,imshow(convexx);
%solidd=solid;
axisarea=major/maxArea12;
I22=imcrop(im2,rect);
glcm = graycomatrix(I22,'NumLevels',256,'Offset',[1 0]);
glcm=glcm(120:220,120:220);
stats = graycoprops(glcm,{'all'});
Con=stats.Contrast;
Homo=stats.Homogeneity;
Ene=stats.Energy;
Cor=stats.Correlation;
-