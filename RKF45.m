function [x_new, z] = RKF45(func, x_old, P, dt, L, R, B, Cd, Cv)
    %{
    INPUTS:
    func: calling the function to evaluate RHS of the ODE (F or
    F_modified).
    x_old: 32*1 vector of solution at previous time step.
    P: 14*1 vector of pressures at previous time step.
    flag: 4*1 vector of +1 & -1 to use the correct sign for the turbolent
    term.
    L: 14*1 vector of inertial parameters.
    R: 14*1 vector of viscous parameters.
    B: 4*1 vector of turbolent parameters (just for heart chambers).
    Cd: 4*1 vector of the dynamic characterstics in each heart chambers.
    Cv: 4*1 vector of the Ejection velocity in each heart chambers.

    OUTPUT:
    z: 32*1 vector of solutions at current time (O(h^5))
    x_new: 32*1 vector of solutions at current time (O(h^4))
    NOTE:
    This function evaluates one step using the Fehlberg method.
    %}
    
    % For RHS evaluation we need the signs
    flag_old = sign([x_old(20:21);x_old(27:28)]);
    
    % Butcher tablau coefficients for the Fehlberg method
    b_O4 = [25/216; 0; 1408/2565; 2197/4104; -1/5];% 4th order accuracy
    b_O5 = [16/135; 0; 6656/12825; 28561/56430; -9/50; 2/55]; % 5th order accuracy      
    a = zeros(5, 5);
    a(1, 1) = 1/4;
    a(1:2, 2) = [3/32; 9/32];
    a(1:3, 3) = [1932/2197; -7200/2197; 7296/2197];
    a(1:4, 4) = [439/216; -8; 3680/513; -845/4104];
    a(1:5, 5) = [-8/27; 2; -3544/2565; 1859/4104; -11/40];
    
    % RKF45 Fehlberg Method
    RKsteps = zeros(32, 6);
    RKsteps(:, 1) = dt* func(x_old, P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 2) = dt* func(x_old + (RKsteps(:, 1) * a(1, 1)), P, flag_old, L, R, B, Cd, Cv); 
    RKsteps(:, 3) = dt* func(x_old + (RKsteps(:, 1:2) * a(1:2, 2)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 4) = dt* func(x_old + (RKsteps(:, 1:3) * a(1:3, 3)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 5) = dt* func(x_old + (RKsteps(:, 1:4) * a(1:4, 4)), P, flag_old, L, R, B, Cd, Cv);
    RKsteps(:, 6) = dt* func(x_old + (RKsteps(:, 1:5) * a(1:5, 5)), P, flag_old, L, R, B, Cd, Cv);

    % Update
    z = x_old + RKsteps(:, 1:6) * b_O5 ; 
    x_new = x_old + RKsteps(:, 1:5) * b_O4;
end
