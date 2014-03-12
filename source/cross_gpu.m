function c = cross_gpu(a,b)
% inputs:   a is a 3 x 1 array
%           b is a 3 x 1 array
% outputs:  c is a 3 x 1 array 

% compute the cross product of 2 vectors: c = a x b
c    = gpuArray.zeros(3, 1);
c(1) = a(2)*b(3) - a(3)*b(2);
c(2) = a(3)*b(1) - a(1)*b(3);
c(3) = a(1)*b(2) - a(2)*b(1);

end

