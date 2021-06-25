close all

c0 = 340*1000;
f = 40000;
omega = 2*pi*f;
k = omega/c0;
lamda = c0/f;

% ŽlŠp
sp_x0 = -80:10:70;
sp_y0 = -80:10:70;
[sp_x,sp_y] = meshgrid(sp_x0,sp_y0);

sp_x = reshape(sp_x,length(sp_y0)^2,1);
sp_y = reshape(sp_y,length(sp_y0)^2,1);
sp_z = ones(length(sp_y0)^2,1);
sig = ones(size(sp_z));
%     
delta_x = 1;
delta_y = 1;
delta_z = 1;

wall_z = 15.1*lamda+sp_z(1);
x = (-50:delta_x:50);
y = (-50:delta_x:50);
z = (wall_z-1:delta_x:wall_z);
[X,Y,Z] = meshgrid(x,y,z);



cp_x = [20,1.4*lamda+20,-20,1.4*lamda-20,-1.4*lamda/2,1.4*lamda/2,-25,1.4*lamda-25];
cp_y = [-30,-30,40,40,5,5,-30,-30];
cp_z = [wall_z,wall_z,wall_z,wall_z,wall_z,wall_z,wall_z,wall_z];

    
%ƒsƒXƒgƒ“‚Ì”¼Œa
a =4.5;

sp_z_w = ones(size(sp_z))*((wall_z-sp_z(1))+wall_z);

for ite = 1:400

    cp_p = zeros(length(cp_x));
    for n = 1:length(cp_x)

        for m = 1:length(sp_x)
             P_im = theory_p_flat_one(k,a,cp_x(n),cp_y(n),cp_z(n),sp_x(m),sp_y(m),sp_z_w(m),-1);
             P0 = theory_p_flat_one(k,a,cp_x(n),cp_y(n),cp_z(n),sp_x(m),sp_y(m),sp_z(m),1);
             cp_p(n) = cp_p(n)+sig(n)*(P0+P_im);
        end
        cp_p(n) = cp_p(n)/length(cp_x)/abs(cp_p(n));
    end
    sig = zeros(size(sp_x));
    for q = 1:length(sp_x)
        for l = 1:length(cp_x)
             P_im = conj(theory_p_flat_one(k,a,cp_x(l),cp_y(l),cp_z(l),sp_x(q),sp_y(q),sp_z_w(q),-1));
             P0 = conj(theory_p_flat_one(k,a,cp_x(l),cp_y(l),cp_z(l),sp_x(q),sp_y(q),sp_z(q),1));
             sig(q) = sig(q) + cp_p(l)*(P_im+P0);
        end
        sig(q) = sig(q)/abs(sig(q));
    end
end

P = zeros(size(X));
for n = 1:length(sp_x)

        P_im = theory_p_flat(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z_w(n),-1);
        %P_im = 0;

        P0 = theory_p_flat(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z(n),1);


        P = P+sig(n)*(P0+P_im);

end



Power = 20*log10(abs(P));

xslice =0;
yslice =[];
zslice = wall_z;

figure(3)
hold on
slice(X,Y,Z,Power,xslice,yslice,zslice)
view(0,90)
ax = gca;
ax.FontSize = 15;
xlabel('x (mm)','FontSize',15);
ylabel('y (mm)','FontSize',15);
zlabel('z (mm)','FontSize',15);
title("Amplitude field")
c = colorbar;
c.Label.String = 'Sound pressure level(dB)';
colormap jet
caxis([-20,10])
shading interp
axis equal
% quiver3(sp_x,sp_y,sp_z,zeros(size(sp_x)),zeros(size(sp_x)),sp_z,0.1)
% quiver3(sp_x,sp_y,sp_z_w,zeros(size(sp_x)),zeros(size(sp_x)),-sp_z,0.1)
hold off
