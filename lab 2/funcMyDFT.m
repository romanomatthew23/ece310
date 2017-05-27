function [X] = funcMyDFT(x, M)
%funcMyDFT computes the DFT of x of length M (DFT length)
%   If M>N then zero pad the rest
%   If M<N truncate to M
%   outputs the DFT to the variable X (capital X)
%determine long or wide vector
[v,q] = size(x);
N = max(v,q);
if(v>=q)
X = zeros(M,1);
end

if(v<=q)
X = zeros(1,M);
end
%the DFT is the summation from 0 to N-1 of x[n]*e^(-j2pikn/N)
%in this case K (the DFT length) is M
%and the function length is N (i.e. x[n]=/0 n = 0 to N-1)

%The first and easiest case is when M=N
if(N==M)
    for K = 1:M        %K goes from 0 to M-1 but matlab has weird indices
        X(K) =  0;
        for n = 1:N %n goes from 0 to N-1 but same thing
            %here notice the shifts by 1 to account for the indices
            X(K) = x(n)*exp(-1*i*2*pi*(K-1)*(n-1)/N) + X(K);  
        end
    end 
end

%truncation
if(N>M)
    for K = 1:M        
        X(K) =  0;
        for n = 1:M %notice we're going 1 to M
            %here notice we're dividing by M
            X(K) = x(n)*exp(-1*i*2*pi*(K-1)*(n-1)/M) + X(K);  
        end
    end 
end

%zero-padding
if(N<M)
    for K = 1:M        %K goes from 0 to M-1 but matlab has weird indices
        X(K) =  0;
        for n = 1:N %no goes 1 to N
            %here we divide by M
            X(K) = x(n)*exp(-1*i*2*pi*(K-1)*(n-1)/M) + X(K);  
        end
    end 
end


end

