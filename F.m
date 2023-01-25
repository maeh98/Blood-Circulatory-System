function y = F(x, P, flag, L, R, B, Cd, Cv)
    %{
    INPUTS:
    x: 32*1 vector of solution at previous time step.
    P: 14*1 vector of pressures at previous time step.
    flag: 4*1 vector of +1 & -1 to use the correct sign for the turbolent
    term.
    L: 14*1 vector of inertial parameters.
    R: 14*1 vector of viscous parameters.
    B: 4*1 vector of turbolent parameters (just for heart chambers).
    Cd: 4*1 vector of the dynamic characterstics in each heart chambers.
    Cv: 4*1 vector of the Ejection velocity in each heart chambers.

    OUTPUT:
    y: 32*1 vector evaluating the right hand side of the ODE at current
    time.
    %}

    %V = x(1:14);
    Q = x(15:28);
    O = x(29:32);

    % Initialise
    y = zeros(32, 1);
    
    % Volumes
    y(1:14) = ([Q(14); Q(1:13)] - Q(1:14));

    % Do fluxes Qe
    y(15:28) = ((1./L(1:14)) .* (-[P(2:14); P(1)] + P(1:14) - R(1:14).*Q(1:14))); %excluding y(28)
    % inside heart nonlinearity
    y(20:21) = y(20:21) - ((1./L(6:7)) .* (B(1:2).*Q(6:7).*abs(Q(6:7))));
    y(27:28) = y(27:28) - ((1./L(13:14)) .* (B(3:4).*Q(13:14).*abs(Q(13:14))));

    % Do valves
    y(29:32) = ((1./Cd(1:4)).*(-O(1:4)+0.5.*(1+tanh(Cv(1:4).*(-[P(7:8); P(14); P(1)]+[P(6:7);P(13:14)])))));
   
end
