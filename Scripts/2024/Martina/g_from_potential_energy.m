function g = g_from_potential_energy(U,q)
    %Function that outputs the gravity vector from the 
    %potential energy
    %
    %input:
    %- U = th potential energy (the sum of the potential energy of the single
    %links)
    %- q = a vertical vector q [q1;q2;q3]

    g=vpa(simplify(jacobian(U,q)))'