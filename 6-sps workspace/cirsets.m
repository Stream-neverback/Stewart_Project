i=0;%i������¼ͬһz������Բ������
%����ͬһz�����µİ뾶
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
%�뾶������ɣ�Rmax��Rmin
%����ͬһz�����µ�Բ��
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
%ͬһz�����µ�Բ���������,CIRset
clear i j;