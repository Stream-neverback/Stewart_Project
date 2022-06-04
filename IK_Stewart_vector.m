function L = IK_Stewart(px,py,pz,roll,pitch,yaw)
%% The position of Bi relative to B
B1=[-169.76;-581.96;49.49];
B2=[130.07;-592.10;49.49];
B3=[588.88;143.98;49.49];
B4=[447.74;408.70;49.49];
B5=[-419.12;438;49.49];
B6=[-577.81;183.41;49.49];
%% The position of Pi relative to P
P1=[-350;-346.41;-46.39];
P2=[350;-346.41;-46.39];
P3=[470;-129.9;-46.39];
P4=[125;476.31;-46.39];
P5=[-125;476.31;-46.39];
P6=[-475;-129.9;-46.39];
%% Calculate Transform matrix
  T = [cos(pitch)*cos(yaw),                                 -cos(pitch)*sin(yaw),                                  sin(pitch),           px;
        sin(roll)*sin(pitch)*cos(yaw)+cos(roll)*sin(yaw),-sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw),-sin(roll)*cos(pitch),py;
       -cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw)+sin(roll)*cos(yaw), cos(roll)*cos(pitch),pz;
        0,                                                      0,                                                      0,                    1 ];
% or
X=[px;py;pz];
RX = [1 0 0;0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
RY = [cos(pitch) 0 sin(pitch);0 1 0;-sin(pitch) 0 cos(pitch)];
RZ = [cos(yaw) -sin(yaw) 0;sin(yaw) cos(yaw) 0; 0 0 1];
R = RZ*RY*RX;
%% Calculate P relative to B
% for i = 1:6
%     P_B(:,i) = T*P_P(:,i);
% end
%% Calculate the length
l1=X+R*P1-B1;
l2=X+R*P2-B2;
l3=X+R*P3-B3;
l4=X+R*P4-B4;
l5=X+R*P5-B5;
l6=X+R*P6-B6;
L=[l1.';l2.';l3.';l4.';l5.';l6.'];
end