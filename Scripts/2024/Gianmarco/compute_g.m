function gq = compute_g(m,r0c,g,q,cartesian_dimension)
    % m -> column or row
    % g -> row
    % p -> each component r0_i as column
    n = length(m);
    U_vec = [];
    g_vec = [];

    for i=1:n
        U_i = simplify(-m(i)*g*r0c(1:cartesian_dimension, i));
        fprintf("U_%d = ",i);
        disp(U_i);        
        U_vec = [U_vec,U_i];
    end

    U = collect(simplify(sum(U_vec)),[sin(q) cos(q)]);
    fprintf("U = sum(U_i) = ");
    disp(U)

    for i=1:n
        g_i = diff(U,q(i));
        g_vec = [g_vec;g_i];
    end
    gq = g_vec;
end

