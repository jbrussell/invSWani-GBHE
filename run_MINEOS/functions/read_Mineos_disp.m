function [ model ] = read_Mineos_disp( data, model, param )
% Read dispersion from Mineos files
FRECHETPATH = param.frechetpath;
CARDID = param.CARDID;
qname_S = [CARDID,'.',param.STYPEID,'.q'];
qname_T = [CARDID,'.',param.TTYPEID,'.q'];
I_R0 = data.rayl.mode_br_ani == 0;
I_R1 = data.rayl.mode_br_ani == 1;
I_L0 = data.love.mode_br_ani == 0;
I_L1 = data.love.mode_br_ani == 1;

%S0
[phV,grV,phVq] = readMINEOS_qfile2([FRECHETPATH,qname_S],data.rayl.periods_ani(I_R0),0);
model.R0.phv = phV';
model.R0.grv = grV';
model.R0.periods = data.rayl.periods_ani(I_R0)';

%S1
[phV,grV,phVq] = readMINEOS_qfile2([FRECHETPATH,qname_S],data.rayl.periods_ani(I_R1),1);
model.R1.phv = phV';
model.R1.grv = grV';
model.R1.periods = data.rayl.periods_ani(I_R1)';

%T0
[phV,grV,phVq] = readMINEOS_qfile2([FRECHETPATH,qname_T],data.love.periods_ani(I_L0),0);
model.L0.phv = phV';
model.L0.grv = grV';
model.L0.periods = data.love.periods_ani(I_L0)';

%T1
[phV,grV,phVq] = readMINEOS_qfile2([FRECHETPATH,qname_T],data.love.periods_ani(I_L1),1);
model.L1.phv = phV';
model.L1.grv = grV';
model.L1.periods = data.love.periods_ani(I_L1)';

end

