function [b,e] = get_standard_deviation(N,S)


global T

x0 = [2,0.4];

options = optimset('TolFun',1e-8,'TolX',1e-8,'MaxIter',1000000,'MaxFunEvals',1000000);

for n=1:N
    shock = normrnd(0,1,[T,S]);
    
    my_obj = @(x)(gen_std(x,shock));
    
    x = fminsearch(my_obj,x0,options);
    
    b(n) = x(1);
    e(n)  = x(2);
    

end

% b_min = min(beta); b_max= max(beta);
% 
% b_step = b_max-b_min/100;
% 
% e_min = min(eta); e_max= max(eeta);



end

