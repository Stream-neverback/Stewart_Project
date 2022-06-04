i=0;%i用来记录同一z坐标下圆的数量
%生成同一z坐标下的半径
for j=1:6
if(roumax(j)^2-(z-w(j))^2>0)
    Rmax(j)=sqrt(roumax(j)^2-(z-w(j))^2);
else
    Rmax(j)=0;
end
if(roumin(j)^2-(z-w(j))^2>0)
    Rmin(j)=sqrt(roumin(j)^2-(z-w(j))^2);
else
    Rmin(j)=0;
end
end
%半径生成完成，Rmax与Rmin
%生成同一z坐标下的圆表
for j=1:6
%CIRset=[];
if((Rmax(j)~=0)&&(Rmin(j)~=0))
    %for j=1:6
        CIRset(:,i+1)=[u(j);v(j);Rmax(j)];
        CIRset(:,i+2)=[u(j);v(j);Rmin(j)];
        %------
        flgset(:,i+1)=1;
        flgset(:,i+2)=-1;
        %------
        i=i+2;
    %end
elseif((Rmax(j)~=0)&&(Rmin(j)==0))
    %for j=1:3
        CIRset(:,i+1)=[u(j);v(j);Rmax(j)];
        %------
        flgset(:,i+1)=1;
        %------
        i=i+1;
    %end
else
end
end
%同一z坐标下的圆表生成完成,CIRset
clear i j;