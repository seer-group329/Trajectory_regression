%�����������˽���ţ���㷨(Beetle Antennae search algorithm, BAS) ��ʵ�ֺ�����ģ�Ϊ�˸����ʵ�ֵ���
function [BestDindex,IterationNum] = BAS(carPose,B5C_x, B5C_y)

n = size(B5C_x);
LimitLeft = 1;
LimitRight = n(2);
IterationPoint = randi([LimitLeft,LimitRight],[1,6]);
Left = min(IterationPoint);
Right = max(IterationPoint);
DisLeft = norm([carPose(1), carPose(2)]-[B5C_x(Left), B5C_y(Left)]);
DisRight = norm([carPose(1), carPose(2)]-[B5C_x(Right), B5C_y(Right)]);
FindOpti = 0;
BestDindex = 0;
IterationNum = 1;
Step = 2;  %����ÿ�ε�����������������

while FindOpti == 0
    %     DisLeft = norm([carPose(1), carPose(2)]-[B5C_x(Left), B5C_y(Left)]);
    %     DisRight = norm([carPose(1), carPose(2)]-[B5C_x(Right), B5C_y(Right)]);
    n1 = size(IterationPoint, 2);
    minDistaceIndex = IterationPoint(1);
    minDistance = norm([carPose(1), carPose(2)]-[B5C_x(IterationPoint(1)), B5C_y(IterationPoint(1))]);
    
    for i = 2:1:n1
        Dis = norm([carPose(1), carPose(2)]-[B5C_x(IterationPoint(i)), B5C_y(IterationPoint(i))]);
        if Dis < minDistance
            minDistance = Dis;
            minDistaceIndex = IterationPoint(i);
        end
    end
    
    if minDistaceIndex == Left
        LimitLeft = max(1, Left - ((Right-Left)- mod(Right-Left, Step))/Step);
    elseif minDistaceIndex == Right
        LimitRight = min(n(2), Right + ((Right-Left)- mod(Right-Left, Step))/Step);
    else
        LimitLeft = Left + ((minDistaceIndex - Left) - mod(minDistaceIndex - Left,Step))/Step;
        LimitRight = Right - ((Right - minDistaceIndex) - mod(Right - minDistaceIndex,Step))/Step;
    end
    
    if Right - Left < 5 || abs(DisLeft - DisRight) < 0.1
        FindOpti = 1;
        BestDindex = minDistaceIndex;
    end
    
    IterationPoint = randi([LimitLeft,LimitRight],[1,6]);
    Left = min(IterationPoint);
    Right = max(IterationPoint);
    DisLeft = norm([carPose(1), carPose(2)]-[B5C_x(Left), B5C_y(Left)]);
    DisRight = norm([carPose(1), carPose(2)]-[B5C_x(Right), B5C_y(Right)]);
    IterationNum = IterationNum + 1;
    
end

end


