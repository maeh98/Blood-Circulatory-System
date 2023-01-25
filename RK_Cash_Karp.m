function [x_new, z] = RK_Cash_Karp(func, x_old, P, dt, L, R, B, Cd, Cv)
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
    z: 32*1 vector of solutions at current time (O(h^5))
    x_new: 32*1 vector of solutions at current time (O(h^4))
    NOTE:
    This function evaluates one step using the Cash_Karp method, which has 
    the same amount of stages as the Fehlberg method.
    %}
    
    % For RHS evaluation we need the signs
    flag_old = sign([x_old(20:21); x_old(27:28)]);
    
    % Butcher tablau coefficients for the Cash-Karp method
    flag_old = sign([x_old(20:21);x_old(27:28)]);
    b_O4 = [37/378; 0; 250/621; 125/594; 0; 512/1771];% 4th order accuracy
    b_O5 = [2825/27648; 0; 18575/48384; 13525/55296; 277/14336; 1/4]; % 5th order accuracy      
    a = zeros(5, 5);
    a(1, 1) = 1/5;
    a(1:2, 2) = [3/40; 9/40];
    a(1:3, 3) = [3/10; -9/10; 6/5];
    a(1:4, 4) = [-11/54; 5/2; -70/27; 35/27];
    a(1:5, 5) = [ 1631/55296; 175/512; 575/13824; 44275/110592; 253/4096];
    
    % Cash_Karp Method
    RKsteps = zeros(32, 6);
    RKsteps(:, 1) = dt* func(x_old, P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 2) = dt* func(x_old + (RKsteps(:, 1) * a(1, 1)), P, flag_old, L, R, B, Cd, Cv); 
    RKsteps(:, 3) = dt* func(x_old + (RKsteps(:, 1:2) * a(1:2, 2)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 4) = dt* func(x_old + (RKsteps(:, 1:3) * a(1:3, 3)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 5) = dt* func(x_old + (RKsteps(:, 1:4) * a(1:4, 4)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 6) = dt* func(x_old + (RKsteps(:, 1:5) * a(1:5, 5)), P, flag_old, L, R, B, Cd, Cv);

    % Update
    z = x_old + RKsteps(:, 1:6) * b_O5 ; 
    x_new = x_old + RKsteps(:, 1:6) * b_O4;
end
