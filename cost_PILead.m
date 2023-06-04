function [error_PILead] = cost_PILead(controller_params)
    load('param.mat','tf_buck','Percent_DVout');

    fz = controller_params(1);
    fc = controller_params(2);
    pm = controller_params(3);

    Desired_wz = 2*pi*fz;
    Gpi_num = [(1/Desired_wz) 1];
    Gpi_denum = [1 0];
    Gpi = tf(Gpi_num,Gpi_denum);
    G1 = Gpi*tf_buck;

    Desired_wgc = 2*pi*fc;
    
    [mag, phase, wout] = bode(G1, {0.1,1e+15});
    mag = squeeze(mag);
    phase = squeeze(phase);
    wout = squeeze(wout);
    
    G_K1 = interp1(wout,20*log10(mag),Desired_wgc); % at Desired_wgc
    PM1 = interp1(wout,phase,Desired_wgc); % at Desired_wgc
    
    K1 = 10^(G_K1/20);
    Kreq = 1/K1;
    PMreq = -180-PM1+pm;
    K = Kreq*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
    Alpa = Desired_wgc*sqrt((1-(sind(PMreq)))/(1+(sind(PMreq))));
    Beta = Desired_wgc*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
    
    G_Lead_num = [K (K*Alpa)];
    G_Lead_denum = [1 Beta];
    G_Lead = tf(G_Lead_num,G_Lead_denum);
    
    Gpi_lead = G_Lead*Gpi*tf_buck;
    Gpi_lead_CloseLoop = feedback(Gpi_lead, 1);

    error_PILead = error_calc(Gpi_lead_CloseLoop, Percent_DVout);
end

