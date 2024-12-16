function Y = regressor_matrix(arg,a_vector)
    %input
    %a_vector= a vector of form [a1;a2;a3..]
    Y=simplify(jacobian(arg,a_vector))