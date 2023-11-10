function [y,r,p] = Quat_to_EulerZYX(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    y = atan2(2*(q0*q3+q1*q2),1-2*((q2^2)+(q3^2)));
    r = asin(2*(q0*q2-q3*q1));
    p = atan2(2*(q0*q1+q2*q3),1-2*((q1^2)+(q2^2)));
end