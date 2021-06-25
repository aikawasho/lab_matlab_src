close all

global c0
global f
global omega
global k
global a
global sp_x0
global sp_y0
global sp_z0
global trap_p

tp = -10;
trap_p =[0,0,tp];
% éläp
sp_x = -40:10:40;
sp_y = -40:10:40;
sp_z = ones(1,length(sp_y));
[sp_x0,sp_y0,sp_z0] = meshgrid(sp_x,sp_y,sp_z);

sp_x0 = reshape(sp_x0,length(sp_y)^3,1);
sp_y0 = reshape(sp_y0,length(sp_y)^3,1);
sp_z0 = reshape(sp_z0,length(sp_y)^3,1);

% %â~
% phi = 0:2*pi/12:2*pi-pi/12;
% sp_x0 = zeros(96,1);
% sp_y0 = zeros(96,1);
% for l = 1:8
%     
%     sp_x0(1+(l-1)*12:l*12) = (40+(l-1)*10)*cos(phi+pi/10/2*(l-1));
%     sp_y0(1+(l-1)*12:l*12) = (40+(l-1)*10)*sin(phi+pi/10/2*(l-1));
% end
% sp_z0 = ones(96,1);
%     
delta_x = 2;
delta_y = 2;
delta_z = 2;
x = (-100:delta_x:100);
y = (-100:delta_x:100);
z = (-100:delta_x:0);
[X,Y,Z] = meshgrid(x,y,z);

c0 = 340*1000;
f = 40000;
omega = 2*pi*f;
k = omega/c0;

%ÉsÉXÉgÉìÇÃîºåa
a =4.5;
P = zeros(size(X));
xx = 0;
p_n = 1;

phi0 = zeros(8,1);
options = optimoptions(@fminunc,'Display','iter','FunvalCheck','on','MaxFunctionEvaluations',3000,'PlotFcn','optimplotfval');
% filename = 'PlotFcns0818_1'+string(tp)+'.fig';

fun = @object_func2;
[phix,fval] = fminunc(fun,phi0,options);
% ÉOÉâÉfÅ[ÉVÉáÉì
% phix = zeros(8,1);
% for m = 1:8
%     
%     phix(m) = -pi/7*(8-m);
% end
    

for n = 1:length(sp_x0)

  
    P0 = theory_p_flat(k,a,X,Y,Z,sp_x0(n),sp_y0(n),sp_z0(n));


    P = P+P0*exp(1j*(pi*xx+phix(p_n)));
    
    if rem(n,6) == 0
        xx = xx+1;
    end
       
    if rem(n,12) == 0
        p_n = p_n+1;
    end
   
end

Power = 20*log10(abs(P));

FontSize = 15;
xslice =0;
yslice =[];
zslice =[-0,-50];

figure(1)
slice(X,Y,Z,Power,xslice,yslice,zslice)
% colormap(jet);
view(90,0)
ax = gca;
ax.FontSize = FontSize;
xlabel('x (mm)','FontSize',FontSize);
ylabel('y (mm)','FontSize',FontSize);
zlabel('z (mm)','FontSize',FontSize);
title("Amplitude field")
c = colorbar;
c.Label.String = 'Sound pressure level(dB)';
shading interp
axis equal
hold on
quiver3(sp_x0,sp_y0,sp_z0,zeros(size(sp_x0)),zeros(size(sp_x0)),-sp_z0,0.5)
hold off

figure(2)
U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
slice(X,Y,Z,U,xslice,yslice,zslice)


view(90,0)
xlabel('x (mm)','FontSize',FontSize);
ylabel('y (mm)','FontSize',FontSize);
zlabel('z (mm)','FontSize',FontSize);

ax = gca;
ax.FontSize = FontSize;
title("Potential field")
shading interp
c = colorbar;
c.Label.String = 'The GorÅfkov potential';

function f = object_func2(phi)

    global omega
    global c0
    global k
    global a
    global sp_x0
    global sp_y0
    global sp_z0
    global trap_p

    
    p_n = 1;
    delta_x = 1;
    delta_y = 1;
    delta_z = 1;

    x = (trap_p(1)-delta_x:delta_x:trap_p(1)+delta_x);
    y = (trap_p(2)-delta_x:delta_x:trap_p(2)+delta_x);
    z = (trap_p(3)-delta_x:delta_x:trap_p(3)+delta_x);
    [X,Y,Z] = meshgrid(x,y,z);

    xx = 0;
    p_n = 1;
    P = zeros(size(X));

    for n = 1:length(sp_x0)


        P0 = theory_p_flat(k,a,X,Y,Z,sp_x0(n),sp_y0(n),sp_z0(1));


        P = P+P0*exp(1j*(pi*xx+phi(p_n)));

        if rem(n,6) == 0
            xx = xx+1;
        end

        if rem(n,12) == 0
            p_n = p_n+1;
        end

    end

    U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
%     Uxx = diff(U,2,1);
%     Uyy = diff(U,2,2);
%     Uzz = diff(U,2,3);
    
    L = 6*del2(U);
%     
    f = abs(P(2,2,2))-L(2,2,2);
%     f = wp*abs(P(2,2,2))-wx*Uxx(1,2,2)-wy*Uyy(2,1,2)-wz*Uzz(2,2,1);
    
    
end
