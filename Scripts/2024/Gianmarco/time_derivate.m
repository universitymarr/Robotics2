function f_dt = time_derivate(f, var, d_var)
    f_dt = 0;

    for i = 1:length(var)
        current_var = var(i);
        d_current_var = d_var(i);
        
        f_dt = f_dt + diff(f, current_var)*d_current_var;
    end
end