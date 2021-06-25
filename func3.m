function [c,ceq] = func2(phi)

    global wall_z
    global omega
    global c0
    global k
    global a
    global im_z
    global sp_x
    global sp_y
    global sp_z
    global tp
    global wx
    global wy
    global wz
    global wp
    global theta_sp_num
    global reflect_on
    global reverse
    global A
    
    c= zeros(3,1);
    ceq = zeros(3,1);
    ph_n = 1;
    delta_x = 1;
    delta_y = 1;
    delta_z = 1;
    
    z_del = 5;
    x = (tp(1)-delta_x:delta_x:tp(1)+delta_x);
    y = (tp(2)-delta_y:delta_y:tp(2)+delta_y);
    z = (tp(3)-z_del-delta_z:delta_z:tp(3)+z_del+delta_z);
    [X,Y,Z] = meshgrid(x,y,z);

    
    P = zeros(size(X));
    for n = 1:length(sp_x)
        xx = 0;
        if reverse  == 1
            if sp_y(n) < 0 
                xx =1;
            end
        end
        
        P_im = 0;
        P0 = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);
        
        if reflect_on ==1
            P_im = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
        end
        
        P = P+phi(ph_n,2)*(P0+P_im)*exp(1j*(phi(ph_n,1)+pi*xx));
        if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
            ph_n = ph_n +1;
            if ph_n == theta_sp_num+1
                ph_n = 1;
            end
        end
    end

    U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
%     Uxx = diff(U,2,1);
%     Uyy = diff(U,2,2);
    Uzz = diff(U,1,1);
    
    L = 6*del2(U);
    L_mean = mean(mean(mean(abs(L))));
%     
%     f = abs(P(2,2,2))-L(2,2,2);
%     c(1) = abs(L(2,2,end-1))-L_mean/2;
    c(2) = abs(L(2,2,2))-L_mean/2;

    %c(2) = abs(Uzz(2,2,length(z)-1))-0.01;
    
    
end