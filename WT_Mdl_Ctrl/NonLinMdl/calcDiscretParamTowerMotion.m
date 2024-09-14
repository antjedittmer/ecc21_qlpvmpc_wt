%% Parameters for continous functions
initWorkspace;
[wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, pitch, Cq, Ct ] = ...
    initModel5MWNREL;

%foreaft motion
b0fa = 1/wecs.mt;
a1 = -wecs.Kt/wecs.mt;
a0 = -wecs.Bt/wecs.mt;

%sidewards motion
b0sw = 3/(2*wecs.H) * 1/wecs.mt;
a1sw = -wecs.Ktsw/wecs.mt;
a0sw = -wecs.Btsw/wecs.mt;

%% Continous and discrete models
ss_sw_cont = ss([0 1 ; a1sw,a0sw], [1; b0sw], [1,0],0);
ss_fa_cont = ss([0 1 ; a1,a0], [1; b0fa], [1,0],0);
ss_sw_disc = c2d(ss_sw_cont,DT);
ss_fa_disc = c2d(ss_fa_cont,DT);

Asw_z = ss_sw_disc.A;
Bsw_z = ss_sw_disc.B;

