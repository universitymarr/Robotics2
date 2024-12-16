function alpha = bound_of_g(g,q)
    %Function that outputs the upper bound alpha
    % on g 
    %
    %input:
    %- g = the gravity vector
    %- q = a vertical vector q [q1;q2;q3]
    %
    %output: the value of alpha

    dgdq=jacobian(g,q');
    m=simplify(dgdq'*dgdq)
    eigenvalues=eig(m)
    display('maximum eigenvalue:')
    max_=simplify(max(eigenvalues))
    alpha=simplify(sqrt(max_))
    

