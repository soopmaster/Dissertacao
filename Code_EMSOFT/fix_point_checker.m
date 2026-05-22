function checker = fix_point_checker(x1,x2,P1,P2)

cut_off = 1e-6;

dim = size(P1,1);
x1 = reshape(x1,dim,1);
x2 = reshape(x2,dim,1);
%%calculate the exteme points of first ellipsoidal (x-x1)^TP1(x-x1)
[U,V] = eig(inv(P1));

checker = true;

% plot_ellipsoid(x1,P1,'r')
% hold on;
% plot_ellipsoid(x2,P2,'b')
for i = 1:dim
    point1 = x1 + sqrt(V(i,i))*U(:,i);
    point2 = x1 - sqrt(V(i,i))*U(:,i);
%     scatter(point1(1),point1(2));
%     scatter(point2(1),point2(2));
    if (point1 - x2)'*P2*(point1 - x2) > 1 + cut_off|| (point2 - x2)'*P2*(point2 - x2) > 1 + cut_off
        checker = false;
        break;
    end
end
