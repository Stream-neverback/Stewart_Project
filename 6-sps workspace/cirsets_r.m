i=0;%i������¼ͬһ������Բ������
%����ͬһת���µİ뾶
for j=1:6
if(roumax(j)^2-[u(j)^2+v(j)^2-(u(j)+k*v(j))^2/(1+k^2)]>0)
    Rmax(j)=sqrt(roumax(j)^2-[u(j)^2+v(j)^2-(u(j)+k*v(j))^2/(1+k^2)]);
else
    Rmax(j)=0;
end
if(roumin(j)^2-[u(j)^2+v(j)^2-(u(j)+k*v(j))^2/(1+k^2)]>0)
    Rmin(j)=sqrt(roumin(j)^2-[u(j)^2+v(j)^2-(u(j)+k*v(j))^2/(1+k^2)]);
else
    Rmin(j)=0;
end
end
%�뾶������ɣ�Rmax��Rmin
%����ͬһz�����µ�Բ��
for j=1:6
%CIRset=[];
if((Rmax(j)~=0)&&(Rmin(j)~=0))
    %for j=1:6
        CIRset(:,i+1)=[(u(j)+k*v(j))/sqrt(1+k^2);w(j);Rmax(j)];
        CIRset(:,i+2)=[(u(j)+k*v(j))/sqrt(1+k^2);w(j);Rmin(j)];
        %------
        flgset(:,i+1)=1;
        flgset(:,i+2)=-1;
        %------
        i=i+2;
    %end
elseif((Rmax(j)~=0)&&(Rmin(j)==0))
    %for j=1:3
        CIRset(:,i+1)=[(u(j)+k*v(j))/sqrt(1+k^2);w(j);Rmax(j)];
        %------
        flgset(:,i+1)=1;
        %------
        i=i+1;
    %end
else
end
end
%ͬһ�����Բ���������,CIRset
clear i j;