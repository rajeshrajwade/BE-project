
function contactarea = contactareafunc (im2,BW2,angle1)

[m n]=size(BW2);


%for i=1:5
se = strel('line',3,angle1);
I2 = imdilate(BW2,se);

I3 = imdilate(I2,se);
% figure,imshow(I3);
%end



for i=1:m
    for j=1:n
        if( im2(i,j)>=240)
            kav(i,j)=255;
        else
            kav(i,j)=0;
        end
        
    end
end
level =0.1;
kav = im2bw(kav,level);
kav = bwmorph(kav,'clean');
% figure,imshow(kav);
title('skull');

mn=[3 3];
ste = strel('rectangle', mn);
kav1 = imdilate(kav,ste);
kav1 = imdilate(kav1,ste);
%kav1 = imdilate(kav,ste);

%figure,imshow(kav1);
%title('kavti open');



for i=1:m
    for j=1:n
        if( I2(i,j)>=1 && kav1(i,j)>=1)
                k2(i,j)=1;
        else
            k2(i,j)=0;
        end
        
    end
end
% figure,imshow(k2);
title('contact surface area');
om=0;
for i=1:m
    for j=1:n
        if(  k2(i,j)>=1)
                om=om+1;
        end
    end
end
 contactarea=om;





