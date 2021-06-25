function f = object_fun_w(w)

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
    
    ph_n = 1;
    delta_x = 1;
    delta_y = 1;
    delta_z = 1;

    x = (tp(1)-delta_x:delta_x:tp(1)+delta_x);
    y = (tp(2)-delta_y:delta_y:tp(2)+delta_y);
    z = (tp(3)-delta_z:delta_z:tp(3)+delta_z);
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
        
        P = P+(P0+P_im)*exp(1j*pi*xx)*(w(1,ph_n)+1j*w(2,ph_n));
        if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
            ph_n = ph_n +1;
            if ph_n == theta_sp_num+1
                ph_n = 1;
            end
        end
    end

    U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
    Uxx = diff(U,2,1);
    Uyy = diff(U,2,2);
    Uzz = diff(U,2,3);
    
%     L = 6*del2(U);
%     
%     f = abs(P(2,2,2))-L(2,2,2);
    f = wp*abs(P(2,2,2))-wx*Uxx(1,2,2)-wy*Uyy(2,1,2)-wz*Uzz(2,2,1);
    
    
end