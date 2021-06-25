x = 0:500;
P = zeros(length(x),1);
c0 = 346*1000;
f = 40000;
omega = 2*pi*f;
k = omega/c0;

%ƒsƒXƒgƒ“‚Ì”¼Œa
a =4.5;
for x = 0:500

    P(x+1) = theory_p_one(k,a,0,0,x,0,0,0,100);
end


plot(20*log10(abs(P)))