function [y] = syseqn(x,yn1,a)
%y[n] = a*y[n-1] + x[n]
%y=syseqn(x,yn1,a) 
%   Matlab is annoying because the indices are starting at 1 and the
%   function is starting at 0.

%find the length of x and make y that length
[m,n] = size(x);
h = max(m,n);
y = zeros(1,h);

%Initial conditions
y(1) = a*yn1 + x(1)

%rest of values
for k=2:h
    y(k) = a*y(k-1) + x(k);
end



end

