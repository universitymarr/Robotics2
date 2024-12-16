function [T, A] = DHMatrix(DHTABLE)
% Function that takes the DHTABLE and returns T 
% 
% parameters: 
%   - DHTABLE: a n-vector of vectors in the order: [alpha a d theta]
% and outputs:
%   - T: the product of all the matrices corresponding to each vector of
%   the DHTABLE 
%  - A: a cell array containing all the matrices corresponding to each vector of the DHTABLE
%
% example
%
    T = eye(4);
    nums = size(DHTABLE);
    
    A = cell(1,nums(1));
    
    for i = 1:nums(1)
        line = DHTABLE(i, :);
        R = [cos(line(4)) -cos(line(1))*sin(line(4)) sin(line(1))*sin(line(4)) line(2)*cos(line(4));
             sin(line(4)) cos(line(1))*cos(line(4)) -sin(line(1))*cos(line(4)) line(2)*sin(line(4));
             0 sin(line(1)) cos(line(1)) line(3);
             0 0 0 1;];
        A{i} = R;
        T = T * R;   
    end

    if isa(T, 'sym')
        T = simplify(T);
    end
end