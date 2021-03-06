clc
clear
close all
carPose = [3,5,12*pi/36];
theta = carPose(3);
T = carPose;  %平移坐标

carPose = carPose - T;
R = [cos(theta), sin(theta); -sin(theta), cos(theta)]; %旋转矩阵
%生成贝塞尔路径作为原始路径
p0=[10,1,0];
p5=[10,10,0];
p0 = p0 - T;
p5 = p5 - T;
p01 = R*[p0(1); p0(2)];
p51 = R*[p5(1); p5(2)];
p0 = [p01(1), p01(2), p0(3)];
p5 = [p51(1), p51(2), p5(3)];
dis = norm([p0(1), p0(2)]-[p5(1), p5(2)]);
L = dis/2;
p1 = p0 + [L*cos(p0(3)),L*sin(p0(3)),0];
p2 = p1;
p4 = p5 - [L*cos(p5(3)),L*sin(p5(3)),0];
p3 = p4;
% p1 = p0;
% p2 = p1;
% p4 = p5;
% p3 = p4;
[B5C_x, B5C_y] = B5_C(p0,p1,p2,p3,p4,p5);
figure(1)
plot(B5C_x,B5C_y,'k');
theta1 = -theta;
% R1 = [cos(theta1), sin(theta1); -sin(theta1), cos(theta1)];
R1=R';
len = size(B5C_x,2);
% B5C = [0;0];
for i=1:1:len
    B5C = R1*[B5C_x(i); B5C_y(i)];
    B5C = B5C + [T(1);T(2)];
    B5C_x(i) = B5C(1);
    B5C_y(i) = B5C(2);
end
figure(2)
plot(B5C_x,B5C_y,'k');