function e_vt = e_v(t, T_cp, T_rp, T_)
    %{     
    INPUTS 
    t: time since start of contraction
    T_cp: duration of contraction
    T_rp: duration of relaxation
    T_: time of full cardiac cycle
    
    OUTPUTS
    e_a returns the elasticity in the ventricles (a real value).
    %}

    t = mod(t, T_);
    if (t >= 0) && (t<= T_cp)
        e_vt = 0.5*(1-cos(pi*t/T_cp));
    elseif (t >= T_cp) && (t<= T_cp + T_rp)
        e_vt = 0.5*(1+cos(pi*(t-T_cp)/T_rp));
    else
        e_vt = 0;
    end

end