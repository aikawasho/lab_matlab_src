function U = poten_cal(P, delta_x, delta_y, delta_z, c_0, omega) %c_0は媒質の音速, fは音波の周波数, rho_pは微小物体の密度, rho_0は媒質の密度, V_pは微小物体の体積
    
    V_p = 4/3*pi*(1.5)^3; % 3mm球
    c_p = 900*1000;%ポリスチレン 
    rho_0 = 1.3/(1000^3); %空気
    %rho_p = 1700; % 発砲ポリスチレンは密度と速度がわからないの金属とした（かなりアバウト）   
    rho_p = 29/(1000^3); % 発砲ポリスチレン

    [p_x, p_y, p_z] = gradient(P,delta_x,delta_y,delta_z);  
%     
    laplacianP = abs(p_x).^2+abs(p_y).^2+abs(p_z).^2;
    

    
    %ポテンシャルの計算

    K_1 = 1/4*V_p*(1/(c_0^2*rho_0) - 1/(c_p^2*rho_p));
    K_2 = 3/4*V_p*( (rho_p-rho_0)/(omega^2*rho_0*(rho_0+2*rho_p)));
    U = 2*K_1*(abs(P).^2) - 2*K_2*laplacianP;
    
    
%     hold on 
    
    %Forceの計算
    % matlabのgradientを使用してみる
%     [F_x, F_y, F_z] =  gradient(U,delta_x,delta_y,delta_z);
%     F_x = -F_x; F_y = -F_y; F_z = -F_z;

end