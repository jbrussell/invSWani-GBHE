function [ model ] = model_init( card, par, iz_cut )
% Initialize model structure

% Index depths for inversion
rad = card.rad;
Idep = find(rad >= rad(end)-par.model_depth_G*1000);
Idep = Idep(1:end-iz_cut);
IdepB = find(rad >= rad(end)-par.model_depth_B*1000);
IdepB = IdepB(1:end-iz_cut);
IdepH = IdepB;
IdepC = find(rad >= rad(end)-par.model_depth_E*1000);
IdepC = IdepC(1:end-iz_cut);

% Get Values at node midpoints
rad = flipud(rad(Idep));
model.G.z = midpts(flipud(card.z(Idep))')';
model.B.z = midpts(flipud(card.z(IdepB))')';
model.H.z = midpts(flipud(card.z(IdepH))')';
model.E.z = midpts(flipud(card.z(IdepC))')';
model.card_midpts.z = midpts(flipud(card.z(Idep))')';
model.card_midpts.Vsv = midpts(flipud(card.vsv(Idep))')';%/1000);
model.card_midpts.Vsh = midpts(flipud(card.vsh(Idep))')';%/1000);
model.card_midpts.Vsh_E = midpts(flipud(card.vsh(IdepC))')';%/1000);
model.card_midpts.Vpv = midpts(flipud(card.vpv(Idep))')';%/1000);
model.card_midpts.Vph = midpts(flipud(card.vph(Idep))')';%/1000);
model.card_midpts.eta = midpts(flipud(card.eta(Idep))')';%/1000);
model.card_midpts.rho = midpts(flipud(card.rho(Idep))')';%/1000);
model.card_midpts.rho_E = midpts(flipud(card.rho(IdepC))')';%/1000);
model.card = card;

% Vertical Transverse Isotropy parameters
model.VTI.L = model.card_midpts.rho .* model.card_midpts.Vsv.^2;
model.VTI.N = model.card_midpts.rho_E .* model.card_midpts.Vsh_E.^2; 
model.VTI.C = model.card_midpts.rho .* model.card_midpts.Vpv.^2;
model.VTI.A = model.card_midpts.rho .* model.card_midpts.Vph.^2;
model.VTI.F = model.card_midpts.rho .* model.card_midpts.eta.* (model.card_midpts.Vph.^2 - 2*model.card_midpts.Vsv.^2);

end

