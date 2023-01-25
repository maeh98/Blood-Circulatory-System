function [x_new, dt, dt_new, adap_step_count, bjoerck_count, crossed] = adaptive(x_old, t_old, P, dt, L, R, B, Cd, Cv, crossed, METHOD_num, METHOD_diff)
    %{
    INPUTS:
    x_old: 32*1 vector of solution at previous time step.
    t_old: starting time.
    P: 14*1 vector of pressures at previous time step.
    dt: time step size.
    L: 14*1 vector of inertial parameters.
    R: 14*1 vector of viscous parameters.
    B: 4*1 vector of turbolent parameters (just for heart chambers).
    Cd: 4*1 vector of the dynamic characterstics in each heart chambers.
    Cv: 4*1 vector of the ejection velocity in each heart chambers.
    crossed: a counter of how many time some Q changes sign.
    METHOD_num: function handle of the chosen numerical method
    METHOD_diff: function handle corresponding to the chosen way to address
    the non-differentiability.

    OUTPUT:
    x_new: 32*1 vector of solutions at current time (the accuracy depends
    on the method used).
    dt: time step size used.
    dt_new: time step to be used in the next step.
    adap_step_counter: a counter that keeps track of how many times we
    do a step (if the error estimate is too big we need to repeat the step
    with smaller step size).
    bjoerck_count: a counter that finds how many iterations are needed to
    find the root (if such a root needs to be found).
    crossed: a counter that is positive if the any heart chamber flux
    changed sign within the last ten iterations.
    
    NOTE:
    The function evaluates one step with the chosen method. If the error 
    estimate is too big we need to repeat the step with smaller step size.
    We compute the future step size, once the tolerance condition is
    fulfilled. The second part of the function finds (if needed) where the
    flux of a heart chamber crosses 0 via Anderson-Bjoerck root-finding.
    %}


    heart_flow_index = [20;21;27;28];
    tol_bjorck = 1.0e-4;       % Error tolerance for bjorck
    tol_adaptive = 1.0e-4;     % Error tolerance for adaptive time step
    complete_step = 0;
    adap_step_count = 0;
    bjoerck_count = 0;
    
    facmax = 1.5;   % Maximal step size increase, usually [1.5,5]
    facmin = 0.3; 
    p = 4;
    if isequal(METHOD_num, @RK_Bogacki_Shampine)
        p = 2;
    end
    
    while (not(complete_step))
        [x_new, z] = METHOD_num(METHOD_diff, x_old, P, dt, L, R, B, Cd, Cv);
        DEN = norm(x_new-z,2);
        dt_new = dt*min(facmax, max(facmin, 0.9*nthroot(tol_adaptive/(DEN),p+1)));
        adap_step_count = adap_step_count + 1;

        if (DEN > tol_adaptive)
            dt = dt_new;        % Repeat step 
        else
            complete_step = 1;  % Step is completed
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%% Anderson-Bjoerck %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perfom only if the RHS of the ODE is not modified (i.e. it has
    % non-differentiable points). 
    
    if isequal(METHOD_diff, @F)
        % Detect sign changes in heart flow rates
        sign_change = sign(x_old(heart_flow_index)) ~= sign(x_new(heart_flow_index));
        % If Qold~0 OR Qnew~0, don't do Bjoerck.
        root = ((abs(x_old(heart_flow_index)) < tol_bjorck) + (abs(x_new(heart_flow_index)) < tol_bjorck)) >= 1;
        [~,~,index_change] = find(sign_change .* not(root) .* heart_flow_index);

        if (~isempty(index_change)) % There is at least 1 change of sign AND both extremes are not~ 0
            dt0 = 0;                % Left point
            x_new0 = x_old;         % Left evaluation
            found = 0;
            dt_c = dt;              % "Middle" point
            x_c = x_new;            % "Middle" evaluation

            while found == 0        % Until we find the root via Anderson_bjoerck
                m = 0.5*ones(length(index_change),1);
                for l = 1:length(index_change)          % There is potentially a root for more than one variable.
                    if (abs(x_new(index_change(l))) > tol_bjorck && x_c(index_change(l)) < x_new(index_change(l)))
                        m(l) = 1 - x_c(index_change(l)) ./ x_new(index_change(l));
                    end
                end
                
                troot_c = (x_new(index_change) .* (t_old + dt0) - m .* x_new0(index_change) .* (t_old+dt)) ./ (x_new(index_change) - m .* x_new0(index_change));
                [dt_c, closest_root] = min(troot_c - t_old);     % In case more than one heart flux has a root, take the closest one.
                [x_c, ~] = METHOD_num(METHOD_diff, x_old, P, dt_c, L, R, B, Cd, Cv);  % New evaluation.
                bjoerck_count = bjoerck_count + 1;
                
                if (abs(x_c(index_change(closest_root))) < tol_bjorck || dt_c <1e-5)   % Avoid excessively small steps.
                    found = 1;          % End root-finding
                    x_new = x_c;
                    dt = dt_c;                 
                elseif (sign(x_c(index_change(closest_root))) == sign(x_new0(index_change(closest_root))))
                    x_new0 = x_c;       % Use as new left point
                    dt0 = dt_c;
                else
                    x_new = x_c;        % Use as new right point
                    dt = dt_c;
                    
                end
            end
        % Can't use adaptive time step if we cross a point of
        % non-differentiability. Use a fixed step size.
        dt_new = 1.3e-4;
        end
    end
        


end