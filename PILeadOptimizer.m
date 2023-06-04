function [optimum_fz, optimum_fc, optimum_pm] = PILeadOptimizer(StepNum)
    load('param.mat', 'L', 'C', 'fsw');

    fun = @cost_PILead;
    wn = sqrt(1/(L*C));

    rand_fz_values = randi([1 round(wn/(8*pi))],1,StepNum);
    rand_fc_values = randi([round(fsw/20) round(fsw/4)],1,StepNum);
    rand_pm_values = randi([55 70],1,StepNum);

    optimum_fval = inf;
    optimum_param_values = [0 0 0];

    wait_bar = waitbar(0,"PI-Lead Compensator Parameters(fzero, fcrossover and Phase Margin) Optimization in Progress");
    for step=1:StepNum
    [x,fval] = simulannealbnd(fun, [rand_fz_values(step), rand_fc_values(step), rand_pm_values(step)], [1, fsw/20, 55], [wn/(8*pi), fsw/4, 70]);
        if optimum_fval > fval
            optimum_fval = fval;
            optimum_param_values = x;
        end
      waitbar(step/StepNum,wait_bar,"PI-Lead Compensator Parameters(fzero, fcrossover and Phase Margin) Optimization in Progress")
    end
  
    optimum_fz = optimum_param_values(1);
    optimum_fc = optimum_param_values(2);
    optimum_pm = optimum_param_values(3);
    close(wait_bar)
end

