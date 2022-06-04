%% Section 0
clear all;
clc;
close all;

xa=[-169.76 130.07 588.88 447.74 -419.12 -577.81];
ya=[-581.96 -592.10 143.98 408.70 438 183.41];
za=49.49*ones(1,6);
xb=[-350 350 470 125 -125 -475];
yb=[-346.41 -346.41 -129.9 476.31 476.31 -129.9];
zb=-46.39*ones(1,6);
roumin=2201.47*ones(1,6);
roumax=2701.47*ones(1,6);

u=xa-xb;%Բ��λ��x
v=ya-yb;%Բ��λ��y
w=za-zb;%����λ��z

zmin=min(w);%�����ռ���С�߶�
zmax=max(w)+max(roumax);%�����ռ����߶�

jd=1e-9;%�����ж��غϵľ���
ratio=10;%��ɢԲ���ĵ���
ratio_z=16;%��ɢ��z��Ľ���
frame=1;

%% Section 1
%��zֵ��ͼ
%����������������������������������������������������������������������������������������������������������������������������������������������
ARCset_local=[];%��ʼ������
s_local=0;%��ʼ����������
AREA=[];

figure(1);

for z=zmin:2:zmax
    %display(['z = ',num2str(z),', initializing...']);
    CIRset=[];%��ʼ��Բ����
    flgset=[];%��ʼ������߽��ʶ����1Ϊ��-1Ϊ��
    Flagset=[];%��ʼ����ʶ�����ڻ���
    %display('Generating circle sets...');
    cirsets;%����Բ���õ�ͬһz�����µ�Բ��CIRset
    ARCset=[];%��ʼ��������
    %display('Generating arc sets...');
    arcsets;%���ɽ����߽绡,�õ�ͬһz�����µĻ���ARCset
    %��ͬһz�����½����߽绡��ARCset�����ܻ���ARCset_local��
    if(~isempty(ARCset))
        s_ztemp=size(ARCset,2);
        ARCset_local(:,s_local+1:s_local+s_ztemp)=[ARCset;z*ones(1,s_ztemp)];
        AREA=[AREA;calarea(ARCset,Flagset),z];
        s_local=s_local+s_ztemp;
        %-----------------------
        %display('Generating point sets...');
        %PTset=[];
        for i=1:s_ztemp
            t_temp=linspace(ARCset(4,i),ARCset(5,i),ratio);
            r_temp=ARCset(3,i);
            x_temp=ARCset(1,i)+r_temp*cos(t_temp);
            y_temp=ARCset(2,i)+r_temp*sin(t_temp);
            %PTset(:,ratio*(i-1)+1:ratio*i)=[x_temp;y_temp;z*ones(1,ratio)];
            plot3(x_temp,y_temp,z*ones(1,ratio));hold on;axis equal;grid on;%axis vis3d;            
        end        
        %text(0,0,num2str(AREA(1)));
        pause(0.1);
        display(['Local generation done. Area for z=',num2str(z),' is ',num2str(calarea(ARCset,Flagset))]);
        %-----------------------
    else
        display(['z=',num2str(z),' not in the workspace']);
    end
    %ͬһz�����½����߽绡��ARCset�������ܻ���ARCset_local�У�����������¼����
end
clear s_ztemp CIRset ARCset z Rmax Rmin Flagset Flagset_temp PTset flgset s_local t_temp;

%% Section 2
%��z��ƽ���ͼ
%����������������������������������������������������������������������������������������������������������������������������������������������
ARCset_local=[];%��ʼ������
s_local=0;%��ʼ����������
AREA_r=[];

for theta_z=-pi/2+pi/ratio_z/2:pi/ratio_z:pi/2-pi/ratio_z/2;
%for theta_z=-pi/2+pi/ratio_z/2
    k=tan(theta_z);
    CIRset=[];%��ʼ��Բ����
    flgset=[];%��ʼ������߽��ʶ����1Ϊ��-1Ϊ��
    Flagset=[];%��ʼ����ʶ�����ڻ���
    cirsets_r;%������ת�����Բ���õ�ͬһ�����Բ��CIRset
    ARCset=[];%��ʼ��������
    %display('Generating arc sets...');
    arcsets_r;%���ɽ����߽绡,�õ�ͬһ����Ļ���ARCset
    %��ͬһ���潻���߽绡��ARCset�����ܻ���ARCset_local��
    k=tan(theta_z);
    if(~isempty(ARCset))
        s_ztemp=size(ARCset,2);
        ARCset_r_local(:,s_local+1:s_local+s_ztemp)=[ARCset(1,:)/sqrt(1+k^2);...
            k*ARCset(1,:)/sqrt(1+k^2);ARCset(2:5,:);theta_z*ones(1,s_ztemp)];
        s_local=s_local+s_ztemp;
        %-----------------------
        %display('Generating point sets...');
        %PTset=[];
        for i=1:s_ztemp
            t_temp=linspace(ARCset(4,i),ARCset(5,i),ratio);
            r_temp=ARCset(3,i);
            xbar_temp=ARCset(1,i)+r_temp*cos(t_temp);
            x_temp=xbar_temp/sqrt(1+k^2);
            y_temp=xbar_temp*k/sqrt(1+k^2);
            z_temp=ARCset(2,i)+r_temp*sin(t_temp);            
            %PTset(:,ratio*(i-1)+1:ratio*i)=[x_temp;y_temp;z*ones(1,ratio)];
            plot3(x_temp,y_temp,z_temp);hold on;axis equal;grid on;%axis vis3d;            
        end        
        %text(0,0,num2str(AREA(1)));
        pause(0.1);
        display('Local generation done.');
        %-----------------------
    else
        display(['theta=',num2str(theta_z),' not in the workspace']);
    end
    %ͬһz�����½����߽绡��ARCset�������ܻ���ARCset_local�У�����������¼����
end
clear s_ztemp CIRset ARCset z Rmax Rmin Flagset Flagset_temp PTset flgset jd s_local t_temp;