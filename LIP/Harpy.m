%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harpy Class Geometry (2D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Harpy

    % any class attributes
    properties
        L  % link lengths in [m]
        B  % base (hip) frame location in x-z plane [m]
    end

    % any class functions
    methods
        % constructor for the class
        function obj = Harpy(B,L)
            obj.B = B;
            obj.L = L;
        end

        % forward kineamtics
        function p = fwd_kin(obj,q)
            % extract link length info
            L1 = obj.L(1);
            L2 = obj.L(2);
            L3 = obj.L(3);

            % joint angles
            q1 = q(1);
            q2 = q(2);

            % compute joint locations
            p_hip = [0; 0];
            p_knee = [L1*cos(q1); -L1*sin(q1)];
            p_ankle = p_knee + [L2*cos(q1 + q2); -L2*sin(q1+q2)];
            p_foot = p_ankle + [L3*cos(q1); -L3*sin(q1)];
            
            p = [p_hip, p_knee, p_ankle,p_foot] + [obj.B(1); obj.B(2)];
        end

        % compute jacobian
        function J = jacobian(obj,q)
            % extract link length info
            L = obj.L;
            L1 = L(1);
            L2 = L(2);
            L3 = L(3);

            % joint angles
            q1 = q(1);
            q2 = q(2);
        
            J(1,1) = -(L1+L3)*sin(q1) - L2*sin(q1+q2);
            J(1,2) = -L2*sin(q1+q2);
            J(2,1) = -(L1+L3)*cos(q1) - L2*cos(q1+q2);
            J(2,2) = -L2*cos(q1+q2);
        end

        % compute inverse kinematics, Newton raphson method
        function q = inv_kin(obj,p,qk)
           
            % IK settings
            err_tol = 0.001;
            max_iter = 500;
        
            % iterate until convergence
            err = Inf;
            iters = 0;
            while (err > err_tol)

                % compute error so far
                xk = obj.fwd_kin(qk);
                e = p - xk(:,end);
        
                % jacobian of old guess
                J = obj.jacobian(qk);
        
                % iterate for a new guess
                qk = qk + inv(J) * e;
        
                err = norm(e);
                iters = iters+1;
                if iters > max_iter
                    break
                end
            end
            
            % return last iterated point
            q = qk;
        end
    end
end