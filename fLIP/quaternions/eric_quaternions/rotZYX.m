function [v1] = rotZYX(v0,yaw,roll,pitch)
    c1 = cos(yaw); c2 = cos(roll); c3 = cos(pitch);
    s1 = sin(yaw); s2 = sin(roll); s3 = sin(pitch);
    
    R = [c1*c2,c1*s2*s3-c3*s1,s1*s3+c1*c3*s2;...
         c2*s1,c1*c3+s1*s2*s3,c3*s1*s2-c1*s3;...
         -s2  ,c2*s3         ,c2*c3         ];
     
    v1 = R*v0;
end