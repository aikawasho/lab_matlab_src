x = -20:0.1:20;
[X,Y,Z] = meshgrid(x,x,x);

myu = [0;0;0];

delta = [200,0,0;0,200,0;0,0,100];

A = zeros(size(X));

B = zeros(length(x));
B2 = ones(length(x),1);
B2((x>=-10)) = 5;
B2((x>=-5)) = 10;
B2((x>5)) = 5;
B2((x>=10)) = 1;


for i = 1:length(X)
    B(i) = 1/sqrt(2*pi*70)*exp(-x(i)^2/2/70);
end
B = B/max(B)*10;
figure(1)
plot(x,B)
hold on
grid on
plot(x,B2)

hold off

% [X2,Y2] = meshgrid(x,x);
% C = zeros(size(X2));
% delta2 = [100,0;0,100];
% myu2 = [0;0];
% 
% for i=1:length(X)
%     for j = 1:length(X)
%         tmp =[x(i);x(j)];
%         C(j,i) = 1/sqrt(2*pi)^(2/2)/sqrt(det(delta2))*exp(-1/2*(tmp-myu2)'*delta2^(-1)*(tmp-myu2));
%     end
% end
% figure(2)
% surf(X2,Y2,C)
% shading interp
% xlabel('x (mm)');
% ylabel('y (mm)');
% for i = 1:length(X)
%     for j = 1:length(X)
%         for k = 1:length(X)
%             tmp = [x(i);x(j);x(k)];
%             A(j,i,k) = 1/sqrt(2*pi)^(4/2)/sqrt(det(delta))*exp(-1/2*(tmp-myu)'*delta^(-1)*(tmp-myu));
%         end
%     end
% end
% 
% figure(3)
% xslice =0;
% yslice =[];
% zslice = 0;
% slice(X,Y,Z,A,xslice,yslice,zslice)
% xlabel('x (mm)');
% ylabel('y (mm)');
% zlabel('z (mm)');
% shading interp
