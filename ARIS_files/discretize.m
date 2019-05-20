function [Ad,Bd,Kd] = discretize(A, B, K, Ts)
    Ad = expm(A*Ts);
    Bd = A\(Ad-eye(size(A)))*B;
    Kd = A\(Ad-eye(size(A)))*K;
end