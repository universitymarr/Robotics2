function is_an_inertia_matrix = is_inertia_matrix(matrix)
    
    % Check if the matrix is symmetric
    if ~isequal(matrix, matrix')
        error('Matrix is not symmetric');
    end
    
    % Check if the matrix is positive definite
    % Check determinant
    if isa(matrix, 'sym')
        warning('The matrix is symbolic. Please check manually that the determinant is not negative.');
        disp("determinant:")
        disp(det(matrix))
    elseif det(matrix) <= 0
        is_an_inertia_matrix = false;
        disp("this matrix is not an inertia matrix since the determinant is negative")
        return;
    end
    
    % Compute eigenvalues
    eigenvalues = eig(matrix);
    
    % Check if all eigenvalues are positive
    if isa(matrix, 'sym')
        warning('The matrix is symbolic. Please check manually that the eigenvalues are all positive.');
        disp("eigenvalues:")
        disp(eigenvalues)
        is_an_inertia_matrix = [];
    elseif all(eigenvalues > 0)
        is_an_inertia_matrix = true;
    else
        is_an_inertia_matrix = false;
        disp("this matrix is not an inertia matrix since some eigenvalues are negative")
    end
    
    if is_an_inertia_matrix == true
        disp("this is an inertia matrix, remember to check that the matrix is not dependent directly on q1!")
    end
end