function [u1,u2,u3] = controlAttitude(qY,qR,qP,wY,wR,wP,q_des,theta)
    kp = 200;
    kd = 4;
    
    E  = [qY;qR;qP];
    E0 = q_des;
    w  = [wP,wR,wY];
    
    a  = rotZYX([0;0;1],E(1),E(2),E(3));
    a0 = rotZYX([0;0;1],E0(1),E0(2),E0(3));
    v1 = rotZYX([0;cos(theta);sin(theta)],E(1),E(2),E(3));
    v2 = rotZYX([-(sqrt(3)/2)*cos(theta);-0.5*cos(theta);sin(theta)],E(1),E(2),E(3));
    v3 = rotZYX([ (sqrt(3)/2)*cos(theta);-0.5*cos(theta);sin(theta)],E(1),E(2),E(3));
    
    q  = eulerZYX_to_Quat(E(1),E(2),E(3));
    q0 = eulerZYX_to_Quat(E0(1),E0(2),E0(3));
    qe = quatMultiply(q0,quatConj(q));
    
    [v,~] = quatAxis(qe);
    Aerr  = qe(2:4);
    
    u1 = kp*norm(Aerr,2)*dot(v1,v) + kd*norm(w,2)*dot(w,v1);
    u2 = kp*norm(Aerr,2)*dot(v2,v) + kd*norm(w,2)*dot(w,v2);
    u3 = kp*norm(Aerr,2)*dot(v3,v) + kd*norm(w,2)*dot(w,v3);
end