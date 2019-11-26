function [ cij_calc ] = cij_calculations( cij )
% JBR - 6/30/17
%
% Calculate cij values
% - G,C,B,H
% - A,C,F,L,N
% - Strength G,C,B
% - Direction G,C,B
%
% JBR 9/13/17 - Removed multiplication by 2 when calculating G/L, B/A, etc
%
% JBR 11/29/17 - Allow all cij's to be calculated at once
%
if size(cij,1) ~= 1
    niter = 1;
    single_cij = 1;
else
    niter = length(cij);
    single_cij = 0;
end

for ic = 1:niter
    
    if single_cij
        c = cij;
        ref = '';
    else
        c = cij(ic).c;
        ref = cij(ic).ref;
    end
    
    % 2theta terms
    bc = 0.5.*(c(1,1)-c(2,2));
    bs = c(1,6) + c(2,6);
    gc = 1./2.*(c(5,5)-c(4,4));
    gs = c(5,4);
    hc = 1./2.*(c(1,3)-c(2,3));
    hs = c(3,6);

    % 4theta terms
    cc = 1./8.*(c(1,1)+c(2,2)) - 1./4.*c(1,2) - 1./2.*c(6,6);
    cs = 1./2.*(c(1,6)-c(2,6));

    % TI terms
    a = 3./8.*(c(1,1)+c(2,2))+1./4.*c(1,2)+1/2*c(6,6);
    capc = c(3,3);
    f = 1./2.*(c(1,3)+c(2,3));
    l = 1./2.*(c(4,4)+c(5,5));
    n = 1./8.*(c(1,1)+c(2,2))-1./4.*c(1,2) + 1./2.*(c(6,6));

    % Normalize
    bc_a = bc/a;
    bs_a = bs/a;
    gc_l = gc/l;
    gs_l = gs/l;
    hc_f = hc/f;
    hs_f = hs/f;
    cc_n = cc/n;
    cs_n = cs/n;


    % Strength
    strength_g = sqrt(gc.^2+gs.^2);
    strength_g = strength_g/l*100;

    strength_c = sqrt(cc.^2+cs.^2);
    strength_c = strength_c/n*100;

    strength_b = sqrt(bc.^2+bs.^2);
    strength_b = strength_b/a*100;

    % Direction
    fastdir_g = atan2d(gs,gc)/2;
    % if fastdir_g < 0
    %     fastdir_g = fastdir_g + 180;
    % end

    fastdir_c = atan2d(cs,cc)/4;
    % if fastdir_c < 0
    %     fastdir_c = fastdir_g + 90;
    % end

    fastdir_b = atan2d(bs,bc)/2;
    % if fastdir_b < 0
    %     fastdir_b = fastdir_b + 180;
    % end

    cij_calc(ic).bc = bc;
    cij_calc(ic).bs = bs;
    cij_calc(ic).gc = gc;
    cij_calc(ic).gs = gs;
    cij_calc(ic).hc = hc;
    cij_calc(ic).hs = hs;
    cij_calc(ic).cc = cc;
    cij_calc(ic).cs = cs;

    cij_calc(ic).a = a;
    cij_calc(ic).capc = capc;
    cij_calc(ic).f = f;
    cij_calc(ic).l = l;
    cij_calc(ic).n = n;

    cij_calc(ic).bc_a = bc_a;
    cij_calc(ic).bs_a = bs_a;
    cij_calc(ic).gc_l = gc_l;
    cij_calc(ic).gs_l = gs_l;
    cij_calc(ic).hc_f = hc_f;
    cij_calc(ic).hs_f = hs_f;
    cij_calc(ic).cc_n = cc_n;
    cij_calc(ic).cs_n = cs_n;

    cij_calc(ic).strength_g = strength_g;
    cij_calc(ic).strength_c = strength_c;
    cij_calc(ic).strength_b = strength_b;

    cij_calc(ic).fastdir_g = fastdir_g;
    cij_calc(ic).fastdir_c = fastdir_c;
    cij_calc(ic).fastdir_b = fastdir_b;
    
    cij_calc(ic).ref = ref;
end

end

