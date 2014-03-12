function [c1, c2, c3] = crossf(a1, a2, a3, b1, b2, b3)

% compute the cross product of 2 vectors: c = a x b
c1 = a2*b3 - a3*b2;
c2 = a3*b1 - a1*b3;
c3 = a1*b2 - a2*b1;

end

