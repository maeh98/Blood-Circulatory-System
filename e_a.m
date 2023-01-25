function e_at = e_a(t, T_cp, T_rp, t_c, t_r, T_)
    %{     
    INPUTS 
    t: time since start of contraction
    T_cp: duration of contraction
    T_rp: duration of relaxation
    t_c: start of contraction time
    t_r: start of relaxation time 
    T_: time of full cardiac cycle
    
    OUTPUTS
    e_a returns the elasticity in the atria (a real value).
    %}
    
    t = mod(t, T_);
    if (t >= 0) && (t<= t_r + T_rp - T_)
        e_at = 0.5*(1 + cos(pi*(t+T_-t_r)/T_rp));
        
    elseif (t > t_r + T_rp - T_) && (t<= t_c)
        e_at = 0;
        
    elseif (t > t_c) && (t<= t_c+T_cp)
        e_at = 0.5*(1-cos(pi*(t-t_c)/T_cp));
        
    elseif (t > t_c+T_cp) && (t<= T_)
        e_at = 0.5*(1+cos(pi*(t-t_r)/T_rp));
    end

end