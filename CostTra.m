function [BestDindex,IterationNum, A, B] = CostTra(carPose,path_x,path_y,dis, minDindex)

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
Rmin = 0.5;
FindOpti = 0;
BestDindex = [];
PathLen = [];
A = [];
B = [];
IterationNum = 1;
abest = 0;
bbest = 0;
for Index = minDindex:1:LenPath-2
    KAmin = 1e10;
    if carPose(3) == pi/2
        a = 0;
    else
        for a = 0.1:0.5:10
            X = [carPose(1),a,0,path_x(Index),path_dx(Index),path_ddx(Index)];
            Y = [carPose(2),a*tan(carPose(3)),0,path_y(Index),path_dy(Index),path_ddy(Index)];
            [x,y] = polynomialConnection(X, Y);
            if max(abs(y - carPose(2))) > 2*dis
                continue;
            end
            Len = size(x, 2);
            t = 0:1/(Len-1):1;
            poly_dx = diff(x)./diff(t);
            poly_dy = diff(y)./diff(t);
            t1 = t(1:end-1);
            poly_ddx = diff(poly_dx)./diff(t1);
            poly_ddy = diff(poly_dy)./diff(t1);
            poly_dy = poly_dy(1:end-1);
            poly_dx = poly_dx(1:end-1);
            k = ((poly_dx.*poly_ddy)-(poly_dy.*poly_ddx))./(sqrt((poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy)));
            klimit = 1/Rmin;
            if max(abs(k)) >= klimit
                continue;
            end
            len = size(t1,2);
            Ka = sum(k)/len;
            if Ka < KAmin
                KAmin = Ka;
                abest = a;
            end
        end
        a = abest;
    end
%     if a == 0
%         for b = 0.1:0.5:10
%             X = [carPose(1),a,0,path_x(Index),path_dx(Index),path_ddx(Index)];
%             Y = [carPose(2),b,0,path_y(Index),path_dy(Index),path_ddy(Index)];
%             [x,y] = polynomialConnection(X, Y);
%             if max(abs(y - carPose(2))) > 2*dis
%                 continue;
%             end
%             Len = size(x, 2);
%             t = 0:1/(Len-1):1;
%             poly_dx = diff(x)./diff(t);
%             poly_dy = diff(y)./diff(t);
%             t1 = t(1:end-1);
%             poly_ddx = diff(poly_dx)./diff(t1);
%             poly_ddy = diff(poly_dy)./diff(t1);
%             poly_dy = poly_dy(1:end-1);
%             poly_dx = poly_dx(1:end-1);
%             k = ((poly_dx.*poly_ddy)-(poly_dy.*poly_ddx))./(sqrt((poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy)));
%             klimit = 1/Rmin;
%             if max(abs(k)) >= klimit
%                 continue;
%             end
%             len = size(t1,2);
%             Ka = sum(k)/len;
%             if Ka < KAmin
%                 KAmin = Ka;
%                 bbest = b;
%             end
%         end
%         b = bbest;
%     else
%         b = a*tan(carPose(3));
%     end
    b = a*tan(carPose(3));        
    X = [carPose(1),a,0,path_x(Index),path_dx(Index),path_ddx(Index)];
    Y = [carPose(2),b,0,path_y(Index),path_dy(Index),path_ddy(Index)];
    IterationNum =+ 1;
    [x,y] = polynomialConnection(X, Y);
    if max(abs(y - carPose(2))) > 2*dis
        continue;
    end
    pathLen = sum(sqrt(diff(x).*diff(x)+diff(y).*diff(y)));
    Len = size(x, 2);
    t = 0:1/(Len-1):1;
    poly_dx = diff(x)./diff(t);
    poly_dy = diff(y)./diff(t);
    t1 = t(1:end-1);
    poly_ddx = diff(poly_dx)./diff(t1);
    poly_ddy = diff(poly_dy)./diff(t1);
    poly_dy = poly_dy(1:end-1);
    poly_dx = poly_dx(1:end-1);
    k = ((poly_dx.*poly_ddy)-(poly_dy.*poly_ddx))./(sqrt((poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy).*(poly_dx.*poly_dx+poly_dy.*poly_dy)));
    klimit = 1/Rmin;
    if max(abs(k)) >= klimit
        continue;
    else
        FindOpti = 1;
        BestDindex(end+1) = Index;
        PathLen(end+1) = pathLen;
        A(end+1) = a;
        B(end+1) = b;
    end
    
    if(FindOpti && size(BestDindex,2)>100)
        break;
    end
end
i = find(PathLen == min(PathLen));
BestDindex = BestDindex(i);
A = A(i);
B = B(i);
end