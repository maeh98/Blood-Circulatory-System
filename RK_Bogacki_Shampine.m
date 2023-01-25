function [x_new, z] = RK_Bogacki_Shampine(func, x_old, P, dt, L, R, B, Cd, Cv)
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
    x_new: 32*1 vector of solutions at current time (O(h^3))
    z: 32*1 vector of solutions at current time (O(h^2))
    NOTE:
    This function evaluates one step using the Bogacki-Shampine method.
    %}
    
    % For RHS evaluation we need the signs
    flag_old = sign([x_old(20:21);x_old(27:28)]);
    
    % Butcher tablau coefficients for the Bogacki-Shampine method
    b_O2 = [2/9; 1/3; 4/9; 0];% 3rd order accuracy
    b_O3 = [7/24; 1/4; 1/3; 1/8]; % 2nd order accuracy      
    a = zeros(3, 3);
    a(1, 1) = 1/2;
    a(1:2, 2) = [0; 3/4];
    a(1:3, 3) = [2/9; 1/3; 4/9];
   
    % Bogacki–Shampine Method
    RKsteps = zeros(32, 4);
    RKsteps(:, 1) = dt* func(x_old, P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 2) = dt* func(x_old + (RKsteps(:, 1) * a(1, 1)), P, flag_old, L, R, B, Cd, Cv); 
    RKsteps(:, 3) = dt* func(x_old + (RKsteps(:, 1:2) * a(1:2, 2)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 4) = dt* func(x_old + (RKsteps(:, 1:3) * a(1:3, 3)), P, flag_old, L, R, B, Cd, Cv);
    
    % Update 
    x_new = x_old + RKsteps(:, 1:4) * b_O3 ; 
    z = x_old + RKsteps(:, 1:4) * b_O2;
end
