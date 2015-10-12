
function  histthres = histogramthresfunc(im2)




[counts,x]=imhist(im2);
%figure,stem(x,counts);
j=1;

for i=20:240
 h2(j)=j;   
 k(j)=counts(i+1)-counts(i);
 j=j+1;
end
%figure,stem(h2,k);
kmin=min(k);
kmax=max(k);

for i=1:220
 if(k(i)<-100)
 k(i)=-100;
 end
end
%figure,stem(h2,k);
% title('normalized');



for i=1:220
 
 k(i)=(k(i)-kmin);

end

%imshow(k);
%I=image(k);
%imwrite(k,image,'jpg');
%I=imread('image.jpg');
%I=rgb2gray(I);
for i=1:220
 for j=1:220
 ff(i,j)=k(j);
 end 
end



%figure,imshow(uint8(ff));
K =  wiener2(ff,[3 3]);
%figure, imshow(uint8(K));
 j=1;

 for j=1:220
 k2(j)=K(50,j);
 end 

for i=1:220
  k2(i)=(k2(i)+kmin);
end
 for i=1:220
  h3(i)=i;
 end
%figure,stem(h3,k2);
%whos k2
%title('k filtered');
i=1;
for j=100:220
 if (k2(j)>(-10)&&k2(j)<(+10))
     
     x1(i)=j;
     i=i+1;
 end 
end
mm=mean(x1)-20;

if(mm>163)
    raj=163;
else
    raj=mm;
end 
 histthres=raj;


