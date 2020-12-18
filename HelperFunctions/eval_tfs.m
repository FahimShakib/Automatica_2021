function [tfs] = eval_tfs(lure, fgrid)

% Define LTI transfer functions
G_yu  = ss(lure.A, lure.B,    lure.C,         0);
G_yw  = ss(lure.A, lure.L,    lure.C,    lure.D);
G_zu  = ss(lure.A, lure.B,    lure.F,    lure.G);
G_zw  = ss(lure.A, lure.L,    lure.F,    lure.H);

% Evaluate LTI transfer functions along relevant frequency grid
tfs.G_yu_frd  = frd(G_yu,fgrid);
tfs.G_yw_frd  = frd(G_yw,fgrid);
tfs.G_zu_frd  = frd(G_zu,fgrid);
tfs.G_zw_frd  = frd(G_zw,fgrid);
end

