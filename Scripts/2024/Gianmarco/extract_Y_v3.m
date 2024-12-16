function Y = extract_Y_v3(Ya, dynamic_symbols)
    % M: inertia matrix
    % C: christoffel matrix
    % G: potential enrgy matrix
    % joint_acceleration: array with joint acceleration
    % this four items are use to build dynamic model

    % dynamic_model: in the form M*ddq + C + G
    % dynamic_coefficient: an array which contain symbol like [ m1+m3 I2
    % ...]
    % dynamic_symbols: how substitute coefficient, so [a1 a2 ... a_n]
    
    Y = sym([]);
    Y_columns = length(dynamic_symbols);
    dynamic_symbols = reshape(dynamic_symbols, 1, Y_columns);

    for i=1:Y_columns
        
        subs_array = zeros(1, Y_columns);
        subs_array(i) = 1;

        Yi = subs(Ya, dynamic_symbols, subs_array);

        Y = [Y, Yi];
    end
    
end


