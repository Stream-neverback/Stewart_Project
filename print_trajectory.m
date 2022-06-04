% x=x.data;
% y=y.data;
% z=z.data;
% 
% A=[x y z];
% B=unique(A,'rows','stable');
% % A = unique(testA,'rows','stable')


% xd=xd.data;
% yd=yd.data;
% % zd=zd.data;
% zd=ones(10001,1);
% C=[xd yd zd];
% D=unique(C,'rows','stable');


x=out.x.data;
y=out.y.data;
z=out.z.data;

A=[x y z];
B=unique(A,'rows','stable');
% A = unique(testA,'rows','stable')
t=out.x.time;
x=out.x.data;
y=out.y.data;
z=out.z.data;

L=out.L.data;

figure(1)
plot3(x*1000,y*1000,z*1000,'-b','Linewidth',1.2)
% title('仿真摆角跟踪曲线（未加扰动）')
xlabel('X(mm)')
ylabel('Y(mm)')
zlabel('Z(mm)')
legend('末端轨迹') 
axis([-3000 3000 -3000 3000 0 6000]);
grid

figure(2)
plot(t,L*1000)
xlabel('T(s)')
ylabel('L(mm)')
legend('驱动杆1','驱动杆2','驱动杆3','驱动杆4','驱动杆5','驱动杆6') 
axis([0 10 -10 10 ]);

% A=[t x.*1000 y.*1000 z.*1000];
B=[t L*1000]