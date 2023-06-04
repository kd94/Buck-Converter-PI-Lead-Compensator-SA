function [error] = error_calc(tf, DVout)
    [step_response_y, step_response_t] = step(tf);
    alpha = 2; % exp(1)^(-2) = 0.1353

    error_ss = abs(1 - step_response_y(end)); % calculate steady state error from step response

    step_info = stepinfo(step_response_y, step_response_t);

    if step_info.Overshoot == 0 % overshoot yoksa
        error_overshoot = 100;
        error_undershoot = 100;
    elseif step_info.Peak < (DVout + 1) % overshoot Dvouttan kucukse
         error_overshoot = 0;
         error_undershoot =  abs(1 - step_info.SettlingMin);
    else % overshoot Dvouttan buyukse
        error_overshoot = abs(step_info.Peak-1)*10;
        error_undershoot = abs(1 - step_info.SettlingMin);
    end

    rise_time = step_info.RiseTime * 10^3;
    settling_time = step_info.SettlingTime * 10^3;

    error = ((1 - exp(1)^(-alpha)) * (error_ss + error_overshoot + error_undershoot)) + (exp(1)^(-alpha) * (settling_time - rise_time));
end

