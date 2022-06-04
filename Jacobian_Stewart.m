function J = Jacobian_Stewart(px,py,pz,roll,pitch,yaw)
% syms px py pz roll pitch yaw;
L = IK_Stewart_vector(px,py,pz,roll,pitch,yaw);
for i = 1:6
    L_vector(:,i) = L(i,:).';
    u_vector(:,i) = (L_vector(:,i))/norm(L_vector(:,i));
end

% The position of Bi relative to B
P1=[-350;-346.41;-46.39];
P2=[350;-346.41;-46.39];
P3=[470;-129.9;-46.39];
P4=[125;476.31;-46.39];
P5=[-125;476.31;-46.39];
P6=[-475;-129.9;-46.39];
RX = [1 0 0;0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
RY = [cos(pitch) 0 sin(pitch);0 1 0;-sin(pitch) 0 cos(pitch)];
RZ = [cos(yaw) -sin(yaw) 0;sin(yaw) cos(yaw) 0; 0 0 1];
R = RZ*RY*RX;
Pio1 = R*P1; Pio2 = R*P2; Pio3 = R*P3; Pio4 = R*P4; Pio5 = R*P5; Pio6 = R*P6;
%% calculate Jocobian
Jacobian(1,:) = [u_vector(:,1).' (cross(Pio1,u_vector(:,1))).'];
Jacobian(2,:) = [u_vector(:,2).' (cross(Pio2,u_vector(:,2))).'];
Jacobian(3,:) = [u_vector(:,3).' (cross(Pio3,u_vector(:,3))).'];
Jacobian(4,:) = [u_vector(:,4).' (cross(Pio4,u_vector(:,4))).'];
Jacobian(5,:) = [u_vector(:,5).' (cross(Pio5,u_vector(:,5))).'];
Jacobian(6,:) = [u_vector(:,6).' (cross(Pio6,u_vector(:,6))).'];
% J = Jacobian;
J = inv(Jacobian);
end