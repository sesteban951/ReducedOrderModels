function [v,theta] = quatAxis(q)
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
    
    v = [0;0;0];
    
    theta = 2*atan2(sqrt(x^2+y^2+z^2),w);
    v(1)  = x/sqrt(x^2+y^2+z^2);
    v(2)  = y/sqrt(x^2+y^2+z^2);
    v(3)  = z/sqrt(x^2+y^2+z^2);
end