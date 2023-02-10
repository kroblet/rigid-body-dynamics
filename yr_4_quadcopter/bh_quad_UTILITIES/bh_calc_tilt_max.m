function tilt_ang_degs = bh_calc_tilt_max( veh_mass_kg, total_thrust_capacity_N )

    weight_force_N = veh_mass_kg * 9.81;

    % total_thrust_capacity_N * sin(comp_ang) = weight_force_N
    comp_ang       = asin( weight_force_N / total_thrust_capacity_N  );

    comp_ang_degs = comp_ang * 180/pi;

    tilt_ang_degs = 90 - comp_ang_degs;
end

