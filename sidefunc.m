
function  side = sidefunc(BW2,S,idx)



[m n]=size(BW2);
cent=S(idx).Centroid;
centx=cent(1,1);
half=round(n/2);
omm=0;omm1=0;

for i=1:m
    for j=1:half
        if(  BW2(i,j)>=1)
                omm=omm+1;
        end
    end
end

for i=1:m
    for j=half:n
        if(  BW2(i,j)>=1)
                omm1=omm1+1;
        end
    end
end

 

if(centx>half||omm1>omm)
    side=1;
elseif(centx<half||omm1<omm)
    side=0 ;
else
    side=2;
        
    end


% if(side>=1)
%     disp('right');
% else
% disp('left');
% end


