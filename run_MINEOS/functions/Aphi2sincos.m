function [ data ] = Aphi2sincos( data, is_RMS )
% Convert strength and direction (A, phi) measurements of anisotropy to sin
% and cosine (Ac, As) equivalents.
%
% dc/c = A*cos(2*(theta - phi))
% dc/c = Ac*cos(2*theta) + As*cos(2*theta)
%
%
% 11/19
% jbrussell.github

err_pct = 1; %0.6; Multiply errors by a constant (for testing)

% Loop over Rayleigh and Love measurements
flds = {'rayl','love'};
for ifld = 1:length(flds)
    fld = (flds{ifld});
    aniso = data.(fld);
    
    % A, Phi
%     c_iso = aniso.c_iso*1000;
    phi2 = aniso.phi2;
    phi4 = aniso.phi4;
    err_phi2 = aniso.err_phi2;
    err_phi4 = aniso.err_phi4;
    A2 = aniso.A2;%.* c_iso;
    A4 = aniso.A4;%.* c_iso;
    % Use RMS error for magnitude fit or error reported by "fit" tool in matlab
    if is_RMS
        err_A2 = aniso.wRMS_2A/2 * err_pct;
        err_A4 = aniso.wRMS_4A/2 * err_pct;
    else    
        err_A2 = aniso.err_2A;%.* c_iso;
        err_A4 = aniso.err_4A;%.* c_iso;
    end
    
    % Convert to Montagner and Nataf form using partial derivatives to convert
    % data uncertainties...
    data.(fld).c2 = A2.*cosd(2*phi2); %c2_save = c2;
    data.(fld).c4 = A4.*cosd(4*phi4); %c4_save = c4;
    data.(fld).s2 = A2.*sind(2*phi2); %s2_save = s2;
    data.(fld).s4 = A4.*sind(4*phi4); %s4_save = s4;
    data.(fld).err_c2 = sqrt( (err_A2.*cosd(2*phi2)).^2 + (-err_phi2*pi/180*2.*A2.*sind(2*phi2)).^2 ); %c2_std_save = c2_std;
    data.(fld).err_c4 = sqrt( (err_A4.*cosd(4*phi4)).^2 + (-err_phi4*pi/180*4.*A4.*sind(4*phi4)).^2 ); %c4_std_save = c4_std;
    data.(fld).err_s2 = sqrt( (err_A2.*sind(2*phi2)).^2 + ( err_phi2*pi/180*2.*A2.*cosd(2*phi2)).^2 ); %s2_std_save = s2_std;
    data.(fld).err_s4 = sqrt( (err_A4.*sind(4*phi4)).^2 + ( err_phi4*pi/180*4.*A4.*cosd(4*phi4)).^2 ); %s4_std_save = s4_std;
end

end

