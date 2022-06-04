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

u=xa-xb;%圆心位置x
v=ya-yb;%圆心位置y
w=za-zb;%球心位置z

zmin=min(w);%工作空间最小高度
zmax=max(w)+max(roumax);%工作空间最大高度

jd=1e-9;%设置判断重合的精度
ratio=10;%离散圆弧的点数
ratio_z=16;%离散过z轴的截面
frame=1;

%% Section 1
%等z值绘图
%———————————————————————————————————————————————————————————————————————
ARCset_local=[];%初始化弧表
s_local=0;%初始化弧表数量
AREA=[];

figure(1);

for z=zmin:2:zmax
    %display(['z = ',num2str(z),', initializing...']);
    CIRset=[];%初始化圆集合
    flgset=[];%初始化内外边界标识集合1为外-1为内
    Flagset=[];%初始化标识，用于弧表
    %display('Generating circle sets...');
    cirsets;%生成圆表，得到同一z坐标下的圆表CIRset
    ARCset=[];%初始化弧集合
    %display('Generating arc sets...');
    arcsets;%生成交集边界弧,得到同一z坐标下的弧表ARCset
    %将同一z坐标下交集边界弧集ARCset放入总弧表ARCset_local中
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
    %同一z坐标下交集边界弧集ARCset被放入总弧表ARCset_local中，弧表数量记录增加
end
clear s_ztemp CIRset ARCset z Rmax Rmin Flagset Flagset_temp PTset flgset s_local t_temp;

%% Section 2
%过z轴平面绘图
%———————————————————————————————————————————————————————————————————————
ARCset_local=[];%初始化弧表
s_local=0;%初始化弧表数量
AREA_r=[];

for theta_z=-pi/2+pi/ratio_z/2:pi/ratio_z:pi/2-pi/ratio_z/2;
%for theta_z=-pi/2+pi/ratio_z/2
    k=tan(theta_z);
    CIRset=[];%初始化圆集合
    flgset=[];%初始化内外边界标识集合1为外-1为内
    Flagset=[];%初始化标识，用于弧表
    cirsets_r;%生成旋转截面的圆表，得到同一截面的圆表CIRset
    ARCset=[];%初始化弧集合
    %display('Generating arc sets...');
    arcsets_r;%生成交集边界弧,得到同一截面的弧表ARCset
    %将同一截面交集边界弧集ARCset放入总弧表ARCset_local中
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
    %同一z坐标下交集边界弧集ARCset被放入总弧表ARCset_local中，弧表数量记录增加
end
clear s_ztemp CIRset ARCset z Rmax Rmin Flagset Flagset_temp PTset flgset jd s_local t_temp;