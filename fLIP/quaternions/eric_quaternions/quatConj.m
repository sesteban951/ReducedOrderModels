function [qc] = quatConj(q)
    qc = q;
    
    qc(2) = -q(2);
    qc(3) = -q(3);
    qc(4) = -q(4);
end