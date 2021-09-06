function [BestDindex,IterationNum, A] = BASCost(carPose,path_x,path_y,dis, minDindex)

LenPath = size(path_x, 2);
t = 0:1/(LenPath-1):1;
path_dx = diff(path_x)./diff(t);
path_dy = diff(path_y)./diff(t);
path_dx(end+1) = path_dx(end);
path_dy(end+1) = path_dy(end);
path_ddx = diff(path_dx)./diff(t);
path_ddy = diff(path_dy)./diff(t);
path_ddx(end+1) = 0;
path_ddy(end+1) = 0;
Rmin = 1;
LimitLeft = minDindex;
LimitRight = LenPath;
IterationPoint = randi([LimitLeft,LimitRight],[1,6]);
Left = min(IterationPoint);
Right = max(IterationPoint);
FindOpti = 0;
BestDindex = 0;
IterationNum = 1;
Step = 3;
A = 0;
MinTotalCostValIndex = 1e10;
MinTotalCostVal = 1e10;
while FindOpti == 0
    
    n1 = size(IterationPoint);
    n1 = n1(2);
    for i = 1:1:n1
        Kmin = 10^10;
        Amin = 0;
        if(carPose(3)== pi/2)
            a = 0;
            b = 1;
        else
            for a = 0.1:1:10
                X = [carPose(1),a,0,path_x(IterationPoint(i)),path_dx(IterationPoint(i)),path_ddx(IterationPoint(i))];
                Y = [carPose(2),a*tan(carPose(3)),0,path_y(IterationPoint(i)),path_dy(IterationPoint(i)),path_ddy(IterationPoint(i))];
                [x,y] = polynomialConnection(X, Y);
                Len = size(x, 2);
                t = 0:1/(Len-1):1;
                poly_dx = diff(x)./diff(t);
                poly_dy = diff(y)./diff(t);
              
%                 poly_dy(end+1) = path_dy(IterationPoint(i));
%                 poly_dx(end+1) = path_dx(IterationPoint(i));
                t1 = t(1:end-1);
                poly_ddx = diff(poly_dx)./diff(t1);
                poly_ddy = diff(poly_dy)./diff(t1);
                poly_dy = poly_dy(1:end-1);
                poly_dx = poly_dx(1:end-1);
%                 poly_ddy(end+1) = path_ddy(IterationPoint(i));
%                 poly_ddx(end+1) = path_ddx(IterationPoint(i));
                
                k = ((poly_dx.*poly_ddy)-(poly_dy.*poly_ddx))./(sqrt((poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy)));
                klimit = 1/Rmin;
                if(max(abs(k)) >= 10)
                    continue;
                end
                if(abs(max(y)) > 2*dis+carPose(2))
                    continue;
                end
                ka = sum(abs(k));
                if ka < Kmin
                    Kmin = ka;
                    Amin = a;
                end
            end
            a = Amin;
            b = a*tan(carPose(3));
        end
        if a == 0
            for b = 0.1:0.5:10
                X = [carPose(1),a,0,path_x(IterationPoint(i)),path_dx(IterationPoint(i)),path_ddx(IterationPoint(i))];
                Y = [carPose(2),b,0,path_y(IterationPoint(i)),path_dy(IterationPoint(i)),path_ddy(IterationPoint(i))];
                [x,y] = polynomialConnection(X, Y);
                Len = size(x, 2);
                t = 0:1/(Len-1):1;
                dx = diff(x)./diff(t);
                dy = diff(y)./diff(t);
                t1 = t(1:end-1);
                ddx = diff(dx)./diff(t1);
                ddy = diff(dy)./diff(t1);
                dy = dy(1:end-1);
                dx = dx(1:end-1);
                k1 = ((dx.*ddy)-(dy.*ddx))./(sqrt((dx.*dx+dy.*dy).*(dx.*dx+dy.*dy).*(dx.*dx+dy.*dy)));
                if(max(abs(k1)) >= 10)
                    continue;
                end
                if(abs(max(y)) > 2*dis+carPose(2))
                    continue;
                end
                ka = sum(abs(k1));
                if ka < Kmin
                    Kmin = ka;
                    Bmin = b;
                end
                b = Bmin;
        end
        X = [carPose(1),a,0,path_x(IterationPoint(i)),path_dx(IterationPoint(i)),path_ddx(IterationPoint(i))];
        Y = [carPose(2),b,0,path_y(IterationPoint(i)),path_dy(IterationPoint(i)),path_ddy(IterationPoint(i))];
        [x,y] = polynomialConnection(X, Y);
        if(max(abs(y))>2*dis+carPose(2))
            continue;
        end
        PathLen = sqrt(diff(x).*diff(x)+diff(y).*diff(y));
        PathLen = sum(PathLen);
        %DisPointRegerssion = norm([carPose(1), carPose(2)]-[path_x(IterationPoint(i)), path_y(IterationPoint(i))]);
        MinDis = norm([carPose(1), carPose(2)]-[path_x(minDindex), path_y(minDindex)]);
        DisCost = PathLen - MinDis;
        Len = size(x, 2);
        t = 0:1/(Len-1):1;
        dx = diff(x)./diff(t);
        dy = diff(y)./diff(t);
        t1 = t(1:end-1);
%         dy(end+1) = path_dy(IterationPoint(i));
%         dx(end+1) = path_dx(IterationPoint(i));
        ddx = diff(dx)./diff(t1);
        ddy = diff(dy)./diff(t1);
%         ddy(end+1) = path_ddy(IterationPoint(i));
%         ddx(end+1) = path_ddx(IterationPoint(i));
        dy = dy(1:end-1);
        dx = dx(1:end-1);
        k1 = ((dx.*ddy)-(dy.*ddx))./(sqrt((dx.*dx+dy.*dy).*(dx.*dx+dy.*dy).*(dx.*dx+dy.*dy)));
        klimit = 1/Rmin;
        ka = sum(abs(k1));
  
        KAverageCost = ka/Len;
        KMaxCost = max(abs(k1));
        if(KMaxCost>20)
            continue;
        end
        if(KMaxCost > 2*klimit)
            w3 = 0;
            Cost = KMaxCost/klimit;
        end
        
        if(KMaxCost < 2*klimit && KMaxCost > klimit)
            w3 = 1;
            Cost = 2;
        end
        
        if(KMaxCost < klimit)
            w3 = 1;
            Cost = 0;
        end
        
        w2 = 3;
        w1 = 2;
        TotalCostVal = w1*DisCost + w2*KAverageCost +w3*KMaxCost + Cost;
        
        if(TotalCostVal < MinTotalCostVal )
            MinTotalCostVal = TotalCostVal;
            MinTotalCostValIndex = IterationPoint(i);
        end
        
    end
    
    if MinTotalCostValIndex == Left
        LimitLeft = max(minDindex, Left - ((Right-Left)- mod(Right-Left, Step))/Step);
    elseif MinTotalCostValIndex == Right
        LimitRight = min(LenPath, Right + ((Right-Left)- mod(Right-Left, Step))/Step);
    elseif MinTotalCostValIndex == 1e10
        continue;
    else
        LimitLeft = Left + ((MinTotalCostValIndex - Left) - mod(MinTotalCostValIndex - Left,Step))/Step;
        LimitRight = Right - ((Right - MinTotalCostValIndex) - mod(Right - MinTotalCostValIndex,Step))/Step;
    end
    
    if Right - Left < 5 
        FindOpti = 1;
        BestDindex = MinTotalCostValIndex;
        A = a;
    end
    
    IterationPoint = randi([min(LimitLeft,LimitRight),max(LimitLeft,LimitRight)],[1,6]);
    Left = min(IterationPoint);
    Right = max(IterationPoint);
    IterationNum = IterationNum + 1;
    
end

end