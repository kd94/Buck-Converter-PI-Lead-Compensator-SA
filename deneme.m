% Design and Performance Analysis of a PWM DC–DC Buck Converter Using PI–Lead Compensator

% Parameters of DC-DC Buck Converter
Vg = 60;
Vo = 48;
Ro = 0.96; %Load resistance
L = 6.2e-6;
rl = 0.0013;
C = 45e-6;
rc = 0.00172; 
rsw = 0.014; %Switch on-resistance
rd = 0.001; %Diode forward resistance
fsw = 200e3; %Switching frequency
Vm = 1; %Sawtooth peak

D = Vo/Vg;
rx = (D*rsw)+((1-D)*rd);
Io = Vo/Ro;
IL = Io;
Gvd_num = [(Ro*(((rd-rsw)*IL)+Vg)*rc*C)/((Ro+rc)*L*C) (Ro*(((rd-rsw)*IL)+Vg))/((Ro+rc)*L*C)];
Gvd_denum = [1 (L+(C*((rc*Ro)+((rx+rl)*(Ro+rc)))))/((Ro+rc)*L*C) (Ro+rx+rl)/((Ro+rc)*L*C)];
Gvd = tf(Gvd_num,Gvd_denum);
Gpwm = 1/Vm;
T1 = Gpwm*Gvd;
figure(1)
margin(T1)
hold on
grid on

figure(2)
bode(T1);
hold on
grid on

Tuncompensated = feedback(T1, 1);
figure(3)
step(Tuncompensated);
hold on 
grid on

%optimum_fval = 0.3095
wn = sqrt(1/(L*C));
wz_max = wn/4;
Desired_fz = 856;
Desired_wz = 2*pi*Desired_fz;
Gpi_num = [(1/Desired_wz) 1];
Gpi_denum = [1 0];
Gpi_step1 = tf(Gpi_num,Gpi_denum);
G1_step1 = Gpi_step1*T1;
figure(1)
margin(G1_step1)
hold on
grid on

figure(2)
bode(G1_step1);
hold on
grid on

G1_CloseLoop_step1 = feedback(G1_step1, 1);
figure(3)
step(G1_CloseLoop_step1);
hold on 
grid on

Desired_PM = 60; %30-45-60-75-90
fgc_max = fsw/4;
fgc_min = fsw/20;
Desired_fgc = 26772;
Desired_wgc = 2*pi*Desired_fgc;

[GM, PM, Wcg, Wcp] = margin(G1_step1); %wcg Magnitude eğrisinde 0dB'ye karşılık gelen frekans (PM'deki), wcp Phase eğrisinde -180 dereceye karşılık gelen frekans (GM'deki) 
[mag, phase, wout] = bode(G1_step1, {0.1,1e+15});
mag = squeeze(mag);
phase = squeeze(phase);
wout = squeeze(wout);

G_K1 = interp1(wout,20*log10(mag),Desired_wgc); % at Desired_wgc
PM1 = interp1(wout,phase,Desired_wgc); % at Desired_wgc

K1 = 10^(G_K1/20);
Kreq = 1/K1;
PMreq = -180-PM1+Desired_PM;
K = Kreq*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
Alpa = Desired_wgc*sqrt((1-(sind(PMreq)))/(1+(sind(PMreq))));
Beta = Desired_wgc*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));

G_Lead_num = [K (K*Alpa)];
G_Lead_denum = [1 Beta];
G_Lead_step1 = tf(G_Lead_num,G_Lead_denum);

Gpi_lead_step1 = G_Lead_step1*Gpi_step1*T1;
figure(1)
margin(Gpi_lead_step1)
hold on
grid on

figure(2)
bode(Gpi_lead_step1);
hold on
grid on

Gpi_lead_CloseLoop_step1 = feedback(Gpi_lead_step1, 1);
figure(3)
[step_response_y, step_response_t]  = step(Gpi_lead_CloseLoop_step1);
step_info = stepinfo(step_response_y, step_response_t );
error_undershoot = abs(1 - step_info.SettlingMin);
step(Gpi_lead_CloseLoop_step1);
hold on 
grid on

%optimum_fval=0.1351
wn = sqrt(1/(L*C));
wz_max = wn/4;
Desired_fz = 1042;
Desired_wz = 2*pi*Desired_fz;
Gpi_num = [(1/Desired_wz) 1];
Gpi_denum = [1 0];
Gpi_step2 = tf(Gpi_num,Gpi_denum);
G1_step2 = Gpi_step2*T1;
figure(1)
margin(G1_step2)
hold on
grid on

figure(2)
bode(G1_step2);
hold on
grid on

G1_CloseLoop_step2 = feedback(G1_step2, 1);
figure(3)
step(G1_CloseLoop_step2);
hold on 
grid on

Desired_PM = 62; %30-45-60-75-90
fgc_max = fsw/4;
fgc_min = fsw/20;
Desired_fgc = 45808;
Desired_wgc = 2*pi*Desired_fgc;

[GM, PM, Wcg, Wcp] = margin(G1_step2); %wcg Magnitude eğrisinde 0dB'ye karşılık gelen frekans (PM'deki), wcp Phase eğrisinde -180 dereceye karşılık gelen frekans (GM'deki) 
[mag, phase, wout] = bode(G1_step2, {0.1,1e+15});
mag = squeeze(mag);
phase = squeeze(phase);
wout = squeeze(wout);

G_K1 = interp1(wout,20*log10(mag),Desired_wgc); % at Desired_wgc
PM1 = interp1(wout,phase,Desired_wgc); % at Desired_wgc

K1 = 10^(G_K1/20);
Kreq = 1/K1;
PMreq = -180-PM1+Desired_PM;
K = Kreq*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
Alpa = Desired_wgc*sqrt((1-(sind(PMreq)))/(1+(sind(PMreq))));
Beta = Desired_wgc*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));

G_Lead_num = [K (K*Alpa)];
G_Lead_denum = [1 Beta];
G_Lead_step2 = tf(G_Lead_num,G_Lead_denum);

Gpi_lead_step2 = G_Lead_step2*Gpi_step2*T1;
figure(1)
margin(Gpi_lead_step2)
hold on
grid on

figure(2)
bode(Gpi_lead_step2);
hold on
grid on

Gpi_lead_CloseLoop_step2 = feedback(Gpi_lead_step2, 1);
figure(3)
[step_response_y, step_response_t]  = step(Gpi_lead_CloseLoop_step2);
step_info = stepinfo(step_response_y, step_response_t );
error_undershoot = abs(1 - step_info.SettlingMin);
step(Gpi_lead_CloseLoop_step2);
hold on 
grid on

%optimum_fval=0.1278
wn = sqrt(1/(L*C));
wz_max = wn/4;
Desired_fz = 2370;
Desired_wz = 2*pi*Desired_fz;
Gpi_num = [(1/Desired_wz) 1];
Gpi_denum = [1 0];
Gpi_step3 = tf(Gpi_num,Gpi_denum);
G1_step3 = Gpi_step3*T1;
figure(1)
margin(G1_step3)
hold on
grid on

figure(2)
bode(G1_step3);
hold on
grid on

G1_CloseLoop_step3 = feedback(G1_step3, 1);
figure(3)
step(G1_CloseLoop_step3);
hold on 
grid on

Desired_PM = 60; %30-45-60-75-90
fgc_max = fsw/4;
fgc_min = fsw/20;
Desired_fgc = 38027;
Desired_wgc = 2*pi*Desired_fgc;

[GM, PM, Wcg, Wcp] = margin(G1_step3); %wcg Magnitude eğrisinde 0dB'ye karşılık gelen frekans (PM'deki), wcp Phase eğrisinde -180 dereceye karşılık gelen frekans (GM'deki) 
[mag, phase, wout] = bode(G1_step3, {0.1,1e+15});
mag = squeeze(mag);
phase = squeeze(phase);
wout = squeeze(wout);

G_K1 = interp1(wout,20*log10(mag),Desired_wgc); % at Desired_wgc
PM1 = interp1(wout,phase,Desired_wgc); % at Desired_wgc

K1 = 10^(G_K1/20);
Kreq = 1/K1;
PMreq = -180-PM1+Desired_PM;
K = Kreq*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
Alpa = Desired_wgc*sqrt((1-(sind(PMreq)))/(1+(sind(PMreq))));
Beta = Desired_wgc*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));

G_Lead_num = [K (K*Alpa)];
G_Lead_denum = [1 Beta];
G_Lead_step3 = tf(G_Lead_num,G_Lead_denum);

Gpi_lead_step3 = G_Lead_step3*Gpi_step3*T1;
figure(1)
margin(Gpi_lead_step3)
hold on
grid on

figure(2)
bode(Gpi_lead_step3);
hold on
grid on

Gpi_lead_CloseLoop_step3 = feedback(Gpi_lead_step3, 1);
figure(3)
[step_response_y, step_response_t]  = step(Gpi_lead_CloseLoop_step3);
step_info = stepinfo(step_response_y, step_response_t );
error_undershoot = abs(1 - step_info.SettlingMin);
step(Gpi_lead_CloseLoop_step3);
hold on 
grid on

%optimum_fval=0.1019
wn = sqrt(1/(L*C));
wz_max = wn/4;
Desired_fz = 2168;
Desired_wz = 2*pi*Desired_fz;
Gpi_num = [(1/Desired_wz) 1];
Gpi_denum = [1 0];
Gpi_step4 = tf(Gpi_num,Gpi_denum);
G1_step4 = Gpi_step4*T1;
figure(1)
margin(G1_step4)
hold on
grid on

figure(2)
bode(G1_step4);
hold on
grid on

G1_CloseLoop_step4 = feedback(G1_step4, 1);
figure(3)
step(G1_CloseLoop_step4);
hold on 
grid on

Desired_PM = 63; %30-45-60-75-90
fgc_max = fsw/4;
fgc_min = fsw/20;
Desired_fgc = 46446;
Desired_wgc = 2*pi*Desired_fgc;

[GM, PM, Wcg, Wcp] = margin(G1_step4); %wcg Magnitude eğrisinde 0dB'ye karşılık gelen frekans (PM'deki), wcp Phase eğrisinde -180 dereceye karşılık gelen frekans (GM'deki) 
[mag, phase, wout] = bode(G1_step4, {0.1,1e+15});
mag = squeeze(mag);
phase = squeeze(phase);
wout = squeeze(wout);

G_K1 = interp1(wout,20*log10(mag),Desired_wgc); % at Desired_wgc
PM1 = interp1(wout,phase,Desired_wgc); % at Desired_wgc

K1 = 10^(G_K1/20);
Kreq = 1/K1;
PMreq = -180-PM1+Desired_PM;
K = Kreq*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
Alpa = Desired_wgc*sqrt((1-(sind(PMreq)))/(1+(sind(PMreq))));
Beta = Desired_wgc*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));

G_Lead_num = [K (K*Alpa)];
G_Lead_denum = [1 Beta];
G_Lead_step4 = tf(G_Lead_num,G_Lead_denum);

Gpi_lead_step4 = G_Lead_step4*Gpi_step4*T1;
figure(1)
margin(Gpi_lead_step4)
hold on
grid on

figure(2)
bode(Gpi_lead_step4);
hold on
grid on

Gpi_lead_CloseLoop_step4 = feedback(Gpi_lead_step4, 1);
figure(3)
[step_response_y, step_response_t]  = step(Gpi_lead_CloseLoop_step4);
step_info = stepinfo(step_response_y, step_response_t );
error_undershoot = abs(1 - step_info.SettlingMin);
step(Gpi_lead_CloseLoop_step4);
hold on 
grid on

%optimum_fval=0.0972
wn = sqrt(1/(L*C));
wz_max = wn/4;
Desired_fz = 2330;
Desired_wz = 2*pi*Desired_fz;
Gpi_num = [(1/Desired_wz) 1];
Gpi_denum = [1 0];
Gpi_step5 = tf(Gpi_num,Gpi_denum);
G1_step5 = Gpi_step5*T1;
figure(1)
margin(G1_step5)
hold on
grid on

figure(2)
bode(G1_step5);
hold on
grid on

G1_CloseLoop_step5 = feedback(G1_step5, 1);
figure(3)
step(G1_CloseLoop_step5);
hold on 
grid on

Desired_PM = 64; %30-45-60-75-90
fgc_max = fsw/4;
fgc_min = fsw/20;
Desired_fgc = 49906;
Desired_wgc = 2*pi*Desired_fgc;

[GM, PM, Wcg, Wcp] = margin(G1_step5); %wcg Magnitude eğrisinde 0dB'ye karşılık gelen frekans (PM'deki), wcp Phase eğrisinde -180 dereceye karşılık gelen frekans (GM'deki) 
[mag, phase, wout] = bode(G1_step5, {0.1,1e+15});
mag = squeeze(mag);
phase = squeeze(phase);
wout = squeeze(wout);

G_K1 = interp1(wout,20*log10(mag),Desired_wgc); % at Desired_wgc
PM1 = interp1(wout,phase,Desired_wgc); % at Desired_wgc

K1 = 10^(G_K1/20);
Kreq = 1/K1;
PMreq = -180-PM1+Desired_PM;
K = Kreq*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));
Alpa = Desired_wgc*sqrt((1-(sind(PMreq)))/(1+(sind(PMreq))));
Beta = Desired_wgc*sqrt((1+(sind(PMreq)))/(1-(sind(PMreq))));

G_Lead_num = [K (K*Alpa)];
G_Lead_denum = [1 Beta];
G_Lead_step5 = tf(G_Lead_num,G_Lead_denum);

Gpi_lead_step5 = G_Lead_step5*Gpi_step5*T1;
figure(1)
margin(Gpi_lead_step5)
hold on
grid on

figure(2)
bode(Gpi_lead_step5);
hold on
grid on

Gpi_lead_CloseLoop_step5 = feedback(Gpi_lead_step5, 1);
figure(3)
[step_response_y, step_response_t]  = step(Gpi_lead_CloseLoop_step5);
step_info = stepinfo(step_response_y, step_response_t );
error_undershoot = abs(1 - step_info.SettlingMin);
step(Gpi_lead_CloseLoop_step5);
hold on 
grid on