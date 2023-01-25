function [x_new, z] = RK_Dormand_Prince(func, x_old, P, dt, L, R, B, Cd, Cv)
    %{
    INPUTS:
    func: calling the function to evaluate RHS of the ODE (F or
    F_modified).
    x_old: 32*1 vector of solution at previous time step.
    P: 14*1 vector of pressures at previous time step.
    dt: time step size.
    L: 14*1 vector of inertial parameters.
    R: 14*1 vector of viscous parameters.
    B: 4*1 vector of turbolent parameters (just for heart chambers).
    Cd: 4*1 vector of the dynamic characterstics in each heart chambers.
    Cv: 4*1 vector of the Ejection velocity in each heart chambers.

    OUTPUT:
    x_new: 32*1 vector of solutions at current time (O(h^5))
    z: 32*1 vector of solutions at current time (O(h^4))
    NOTE:
    This function evaluates one step using the Dormand-Prince method (one 
    more intermediate step than the Fehlberg method).
    %}
    
    % For RHS evaluation we need the signs
    flag_old = sign([x_old(20:21); x_old(27:28)]);
    
    % Butcher tablau coefficients for the Dormand-Prince method
    b_O4 = [35/384; 0; 500/1113; 125/192; -2187/6784; 11/84; 0];% 4th order accuracy
    b_O5 = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40]; % 5th order accuracy      
    a = zeros(6, 6);
    a(1, 1) = 1/5;
    a(1:2, 2) = [3/40; 9/40];
    a(1:3, 3) = [44/45; -56/15; 32/9];
    a(1:4, 4) = [19372/6561; -25360/2187; 64448/6561; -212/729];
    a(1:5, 5) = [9017/3168; -355/33; 46732/5247; 49/176; -5103/18656];
    a(1:6, 6) = [35/384; 0; 500/1113; 125/192; -2187/6784; 11/84];
    
    % Domand_Prince Method
    RKsteps = zeros(32, 7);
    RKsteps(:, 1) = dt* func(x_old, P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 2) = dt* func(x_old + (RKsteps(:, 1) * a(1, 1)), P, flag_old, L, R, B, Cd, Cv); 
    RKsteps(:, 3) = dt* func(x_old + (RKsteps(:, 1:2) * a(1:2, 2)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 4) = dt* func(x_old + (RKsteps(:, 1:3) * a(1:3, 3)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 5) = dt* func(x_old + (RKsteps(:, 1:4) * a(1:4, 4)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 6) = dt* func(x_old + (RKsteps(:, 1:5) * a(1:5, 5)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 7) = dt* func(x_old + (RKsteps(:, 1:6) * a(1:6, 5)), P, flag_old, L, R, B, Cd, Cv);
    
    % Update
    x_new = x_old + RKsteps(:, 1:7) * b_O5 ; 
    z = x_old + RKsteps(:, 1:7) * b_O4;
end
