k=size(CIRset,2);%kΪԲ����
s_temp=0;%ͬһ�����Բ������
for i=1:k
    %��Բi������Բ�Ľ���
    s=0;
    Pset_temp=[];
    for j=1:i-1
        Pset0=cirints(CIRset(1:2,i),CIRset(3,i),CIRset(1:2,j),CIRset(3,j));
        if(~(isnan(Pset0)))
            Pset0=Pset0(:,Pset0(2,:)>0);
            l=size(Pset0,2);
            Pset_temp(:,s+1:s+l)=Pset0;
            s=s+l;            
        end
    end
    for j=i+1:k
        Pset0=cirints(CIRset(1:2,i),CIRset(3,i),CIRset(1:2,j),CIRset(3,j));
        if(~(isnan(Pset0)))
            Pset0=Pset0(:,Pset0(2,:)>0);
            l=size(Pset0,2);
            Pset_temp(:,s+1:s+l)=Pset0;
            s=s+l;            
        end
    end
    Pset_temp=Pset_temp(:,Pset_temp(2,:)>0);
    %���Բi������Բ�Ľ��㣬���㼯��ΪPset_temp,�������Ϊs
    %���ݽ�����ɢԲi
    if(isempty(Pset_temp))
        ARCset_temp=[CIRset(1:3,i);0;pi];         
        Flagset_temp=flgset(:,i);
        s=1;
    else
        angs=atan2(Pset_temp(2,:)-CIRset(2,i)*ones(1,s),Pset_temp(1,:)-CIRset(1,i)*ones(1,s));
        angs=sort(angs);
        ARCset_temp=[CIRset(1,i)*ones(1,s);CIRset(2,i)*ones(1,s);CIRset(3,i)*ones(1,s);angs;angs(2:s),pi];       
        Flagset_temp=flgset(:,i)*ones(1,s);
    end
    %Բi��ɢΪԲ����Բi��Բ������ΪARCset_temp,Բ������Ϊs
    %Բi��Բ�������жϣ�����Բ�������߽�Ļ�
    ARCset_utemp=ARCset_temp;
    %Flagset_utemp=Flagset_temp;
    for j=1:k
        for j_temp=1:s
            ta_temp=(ARCset_temp(4,j_temp)+ARCset_temp(5,j_temp))/2;
            r_temp=ARCset_temp(3,j_temp);
            x_temp=ARCset_temp(1,j_temp)+r_temp*cos(ta_temp);
            y_temp=ARCset_temp(2,j_temp)+r_temp*sin(ta_temp);
            dlta_temp=(x_temp-CIRset(1,j))^2+(y_temp-CIRset(2,j))^2-(CIRset(3,j))^2;
            if(mod(j,2)==1)
                if(~((dlta_temp<0)||(abs(dlta_temp)<jd)))
                    ARCset_utemp(:,j_temp)=nan;
                    Flagset_temp(:,j_temp)=nan;
                end
            else
            
                if(~((dlta_temp>0)||(abs(dlta_temp)<jd)))
                    ARCset_utemp(:,j_temp)=nan;
                    Flagset_temp(:,j_temp)=nan;
                end
            end
        end
    end
    ARCset_utemp(:,any(isnan(ARCset_utemp)))=[ ];
    Flagset_temp=Flagset_temp(~isnan(Flagset_temp));
    s=size(ARCset_utemp,2);
    %����������Բi�ϵĽ����߽绡���ϣ�����ΪARCset_utemp��������Ϊs
    %����������Բ������ͬһz�����µĻ�����
    ARCset(:,s_temp+1:s_temp+s)=ARCset_utemp;
    Flagset(:,s_temp+1:s_temp+s)=Flagset_temp;
    s_temp=s_temp+s;
    %�õ�ͬһz�����µ�Բ�����߽绡����,��Բi�ϵ�ARCset�Ӽ�
end
%�õ�ͬһz�����µ�Բ�����߽绡����ARCset
clear s i j k ARCset_temp ARCset_utemp Pset_temp j_temp ta_temp r_temp x_temp y_temp dlta_temp Pset0 l angs s_temp;
