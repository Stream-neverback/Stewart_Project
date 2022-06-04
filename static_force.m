clear;

% input factor of x, y, z, yaw, pitch, roll 
x = 1;
y = 1;
z = 1;
yaw = 0;
pitch = 0;
roll = 0;

% input force of each rod;
%F = [10000,10000,10000,10000,10000,10000];
load_temp = 100000;
F = [0;0;load_temp;0;0;0];


% the vector of base points in coordinate bottom
 a1_0 = [0.15, 7/20*sqrt(3), 0]';
 a2_0 = [0.6, -sqrt(3)/10, 0]';
 a3_0 = [0.45, -sqrt(3)/4, 0]';
 a4_0 = [-0.45, -sqrt(3)/4, 0]';
 a5_0 = [-0.6, -sqrt(3)/10, 0]';
 a6_0 = [-0.15, 7/20*sqrt(3), 0]';
 
 % the vector of top points in coordinate of upper flat
 b1_1 = [0.35, sqrt(3)/5, 0]';
 b2_1 = [0.475, 3*sqrt(3)/40, 0]';
 b3_1 = [0.125, -11/40*sqrt(3), 0]';
 b4_1 = [-0.125, -11/40*sqrt(3), 0]';
 b5_1 = [-0.475, 3*sqrt(3)/40, 0]';
 b6_1 = [-0.35, sqrt(3)/5, 0]';
 
% or
P = [x; y; z];
RX = [1 0 0;0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
RY = [cos(pitch) 0 sin(pitch);0 1 0;-sin(pitch) 0 cos(pitch)];
RZ = [cos(yaw) -sin(yaw) 0;sin(yaw) cos(yaw) 0; 0 0 1];
R = RZ*RY*RX;

% calculate ui
u1 = (P + R*b1_1 - a1_0)/norm(P + R*b1_1 - a1_0);
u2 = (P + R*b2_1 - a2_0)/norm(P + R*b2_1 - a2_0);
u3 = (P + R*b3_1 - a3_0)/norm(P + R*b3_1 - a3_0);
u4 = (P + R*b4_1 - a4_0)/norm(P + R*b4_1 - a4_0);
u5 = (P + R*b5_1 - a5_0)/norm(P + R*b5_1 - a5_0);
u6 = (P + R*b6_1 - a6_0)/norm(P + R*b6_1 - a6_0);

% calculate each fi

hat_u1 = [u1; cross((R*b1_1),u1)];
hat_u2 = [u2; cross((R*b2_1),u2)];
hat_u3 = [u3; cross((R*b3_1),u3)];
hat_u4 = [u4; cross((R*b4_1),u4)];
hat_u5 = [u5; cross((R*b5_1),u5)];
hat_u6 = [u6; cross((R*b6_1),u6)];

J = [hat_u1, hat_u2, hat_u3, hat_u4, hat_u5, hat_u6];

F_lod = F\J;

index=1;

for i = 2:6
    if F(i)>F(index)
        index = i;
    end
end


ratio = 10000/F_lod(index);
F_lod = F_lod*ratio;
F = F*ratio;


F_load_max = F(3);


