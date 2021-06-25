function P = theory_p(k,a,X,Y,Z,sp_x,sp_y,sp_z,f_p)
% �ʑ����[��, �U��1�̕��f�����g�̉������v�Z
% k �g��
% a �X�s�[�J���a
% XYZ ����O���b�h(3�����e���\��)
% sp_x,y,x �X�s�[�J�ʒu
%f_p �A���C�œ_
    Dim = sqrt((X-sp_x).^2+(Y-sp_y).^2+(Z-sp_z).^2);
    ip1 =(X-sp_x)*(-sp_x)+(Y-sp_y)*(-sp_y)+(Z-sp_z)*(f_p-sp_z);
    ip2 = Dim*sqrt((0-sp_x).^2+(0-sp_y).^2+(f_p-sp_z).^2);
    theta = acos(ip1./ip2);
    
    dire = bessel_func2(k,a,theta);
    
    P =dire./Dim.*exp(1j*k*Dim);
    %inf��Ԗڂɑ傫���l
    P_m2 = maxk(max(max(abs(P))),2);
    P(P==inf) = P_m2(2);
end