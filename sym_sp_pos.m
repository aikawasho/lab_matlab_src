function [x,y,z] = sym_sp_pos(a,R)

 %a ピストンの直径
% R アレイ半径

    %縦に何こスピーカが並ぶのか
    theta_sp_num = round(2*pi*R/4/a);
    phi_sp_nums = zeros(theta_sp_num,1);
    num =1;
    %横に何個スピーカが並ぶのか
    for theta = 0:pi/2/theta_sp_num:pi/2
        R_tmp = R*sin(theta);
        circ = 2*pi*R_tmp;
        phi_sp_nums(num) = round(((circ)/2-a)/a);
        
        if phi_sp_nums(num) <0
            phi_sp_nums(num) = 0;
        end
            
        num =num+1;
    end
    sp_poss0 = zeros(2*sum(phi_sp_nums),2);

    tmp= 0 ;
    num = 1;
    for theta = 0:pi/2/theta_sp_num:pi/2

        sp_poss0(1+tmp:tmp+phi_sp_nums(num),1) = theta;

        sp_phi0 = pi/phi_sp_nums(num)/2:pi/phi_sp_nums(num):pi-pi/phi_sp_nums(num)/2;

        sp_poss0(1+tmp:tmp+phi_sp_nums(num),2) = sp_phi0;
        tmp = tmp+phi_sp_nums(num);
        sp_poss0(1+tmp:tmp+phi_sp_nums(num),2) = -sp_phi0;
        sp_poss0(1+tmp:tmp+phi_sp_nums(num),1) = theta;
        tmp = tmp+phi_sp_nums(num);
        num = num +1;
    end


    x = R*sin(sp_poss0(:,1)).*cos(sp_poss0(:,2));
    y = R*sin(sp_poss0(:,1)).*sin(sp_poss0(:,2));
    z = R*cos(sp_poss0(:,1));
    
%     figure(2)
    quiver3(x,y,z,-x,-y,-z,0.5)
end

    