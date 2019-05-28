function l = tf_to_latex(sys)
%tf_to_latex Creates latex output of transfer function or state space
%system. Works for SISO and SIMO systems
[num,den] = tfdata(sys);
syms s;
ny = length(sys);
t_sym = ones(ny,1) .* s;
l = strings;
for i=1:ny
    t_sym(i) = poly2sym(num{i},s)/poly2sym(den{i},s);
    l(i) = latex(vpa(t_sym(i),3));
end


