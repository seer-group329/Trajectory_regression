clc
clear
close all
carPose = [8,4,24*pi/72];
p0=[8,1,0];
p5=[8,10,0];
theta = carPose(3);
R = [cos(theta), sin(theta); -sin(theta), cos(theta)]; %��ת����
carPose1 = R*[carPose(1); carPose(2)];
carPose = [carPose1(1), carPose1(2), 0];

%���ɱ�����·����Ϊԭʼ·��
dis = norm([p0(1), p0(2)]-[p5(1), p5(2)]);
L = dis/2;
p1 = p0 + [L*cos(p0(3)),L*sin(p0(3)),0];
p2 = p1;
p4 = p5 - [L*cos(p5(3)),L*sin(p5(3)),0];
p3 = p4;

%����ֱ��
% p1 = p0;
% 
% p2 = p1;
% p4 = p5;
% p3 = p4;
[B5C_x, B5C_y] = B5_C(p0,p1,p2,p3,p4,p5);

% ��·������ת�任
len = size(B5C_x, 2);
for i=1:1:len
    B5C = R*[B5C_x(i); B5C_y(i)];
    B5C_x(i) = B5C(1);
    B5C_y(i) = B5C(2);
end
figure(1)
plot(B5C_x, B5C_y, 'Color', 'r', 'LineWidth', 3);
axis equal
hold on
plot(carPose1(1), carPose1(2), 's','color','k','LineWidth', 2,'markersize',10);
text(carPose1(1), carPose1(2),'carPose','color','b');
LenPath = size(B5C_x, 2);
t = 0:1/(LenPath-1):1;
B5C_dx = diff(B5C_x)./diff(t);
B5C_dx(end+1) = B5C_dx(end);
B5C_dy = diff(B5C_y)./diff(t);
B5C_dy(end+1) = B5C_dy(end);
B5C_ddx = diff(B5C_dx)./diff(t);
B5C_ddy = diff(B5C_dy)./diff(t);
B5C_ddx(end+1) = 0;
B5C_ddy(end+1) = 0;


%�ҵ�С����ԭʼ·������ľ����������Ӷ�ȷ���ع�·���ĵ�����С��ǰ�����򱣳�һ��
[minDindex,IterationNum] = BAS(carPose,B5C_x, B5C_y);


%�ҵ���ѻع��
% [bestIndex, IterationNum1,a]= BASCost(carPose,B5C_x, B5C_y,dis,minDindex);
[bestIndex, IterationNum1,a, b]= CostTra(carPose,B5C_x, B5C_y,dis,minDindex);

%�������Żع�·��
X = [carPose(1),a,0,B5C_x(bestIndex),B5C_dx(bestIndex),B5C_ddx(bestIndex)];
Y = [carPose(2),b,0,B5C_y(bestIndex),B5C_dy(bestIndex),B5C_ddy(bestIndex)];
[xbest,ybest] = polynomialConnection(X, Y);
len = size(B5C_x,2);
% B5C = [0;0];
R1=R';
for i=1:1:len
    B5C = R1*[B5C_x(i); B5C_y(i)];
    B5C_x(i) = B5C(1);
    B5C_y(i) = B5C(2);
end
Len = size(xbest, 2);
for i=1:1:Len
    B5C = R1*[xbest(i); ybest(i)];
    xbest(i) = B5C(1);
    ybest(i) = B5C(2);
end
t = 0:1/(Len-1):1;
figure(2)
plot(xbest, ybest , 'Color', 'b', 'LineWidth', 3);
xlabel('X');
ylabel('Y');
hold on
plot(B5C_x, B5C_y, 'Color', 'r', 'LineWidth', 3);
hold on
plot(xbest(1),ybest(1), 's','color','k','LineWidth', 2,'markersize',10);
text(B5C_x(1),B5C_y(1)+0.5,'Beziser','color','b');
text(xbest(1),ybest(1)-0.4,['(',num2str(xbest(1)),',',num2str(ybest(1)),')'],'color','b')
text(carPose(1), carPose(2)-0.7,'carPose','color','b');
% text(xbest(end/2),ybest(end/2),['(',num2str(xbest(end/2)),',',num2str(ybest(end/2)),')'],'color','b')
text(xbest(50),ybest(50)+.4,'ReturnPath','color','b');
text(xbest(end),ybest(end),['(',num2str(xbest(end)),',',num2str(ybest(end)),')'],'color','b')
text(xbest(end),ybest(end)-.5,'ReturnPoint','color','b');
hold on
plot(B5C_x(minDindex), B5C_y(minDindex), 's','color','k','LineWidth', 2,'markersize',10);
text(B5C_x(minDindex)+.2, B5C_y(minDindex),'minDistancePoint','color','b');
title('�켣�ع�')
figure(3)
subplot(2,2,1)
plot(xbest, ybest , 'Color', 'b', 'LineWidth', 3);
hold on
plot(B5C_x, B5C_y, 'Color', 'r', 'LineWidth', 3);
text(B5C_x(1),B5C_y(1)+0.5,'Beziser','color','b');
text(xbest(1),ybest(1)-0.4,['(',num2str(xbest(1)),',',num2str(ybest(1)),')'],'color','b')
text(carPose(1), carPose(2)-1,'carPose','color','b');
text(xbest(end),ybest(end),['(',num2str(xbest(end)),',',num2str(ybest(end)),')'],'color','b')
text(xbest(end),ybest(end)-1,'ReturnPoint','color','b');
title('�켣�ع�')
% xlim([0,10]);

subplot(2,2,2)
dx = diff(xbest)./diff(t);
dy = diff(ybest)./diff(t);
% dy(end+1) = B5C_dy(bestIndex);
% dx(end+1) = B5C_dx(bestIndex);
t1 = t(1:end-1);
plot(t1, dx, 'color', 'b', 'LineWidth', 2);
hold on
plot(t1, dy, 'color', 'b', 'LineWidth', 2);
title('һ�׵���')
% xlim([0,10]);

subplot(2,2,3)
ddx = diff(dx)./diff(t1);
ddy = diff(dy)./diff(t1);
% ddy(end+1) = B5C_ddy(bestIndex);
% ddx(end+1) = B5C_ddx(bestIndex);
t2 = t1(1:end-1);
plot(t2, ddx, 'color', 'b', 'LineWidth', 2);
hold on
plot(t2, ddy, 'color', 'b', 'LineWidth', 2);
title('���׵���')
% xlim([0,10]);

subplot(2,2,4)
lenk = size(t2,2);
k = [];
for i = 1:1:lenk
    k(end+1) = ((dx(i)*ddy(i))-(dy(i).*ddx(i)))./(sqrt((dx(i).*dx(i)+dy(i).*dy(i))).*((dx(i).*dx(i)+dy(i).*dy(i))).*((dx(i).*dx(i)+dy(i).*dy(i))));
end
plot(t2,k)
text(xbest(1),k(1),['(',num2str(xbest(1)),',',num2str(k(1)),')'],'color','b')
text(xbest(end),k(end),['(',num2str(xbest(end)),',',num2str(k(end)),')'],'color','b')
title('����')
% xlim([0,10]);
