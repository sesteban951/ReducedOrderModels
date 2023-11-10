function [Q] = eulerZYX_to_Quat(yaw,roll,pitch)
    q1 = [cos(yaw/2);0;0;sin(yaw/2)];
    q2 = [cos(roll/2);0;sin(roll/2);0];
    q3 = [cos(pitch/2);sin(pitch/2);0;0];
    
    Q = quatMultiply(quatMultiply(q1,q2),q3);
end