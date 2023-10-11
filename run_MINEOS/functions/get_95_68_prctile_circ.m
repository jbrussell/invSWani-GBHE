function [ x_med, x_u95, x_l95, x_u68, x_l68 ] = get_95_68_prctile_circ( x, fld )
% [ x_med, x_u95, x_l95, x_u68, x_l68 ] = get_95_68_prctile_circ( x, fld )
%
% Get the upper and low 95th and 68th percentiles and median;
%
% This version calculates the circular median and 1 and 2 sigma
%
% x input in degrees
%

if contains(fld,'E') || contains(fld,'4')
    th = 4;
else
    th = 2;
end

x_med = circ_median(th*x*pi/180) * 180/pi / th;
% x_u95 = (x_med + circ_confmean(th*x*pi/180,1-0.95)'* 180/pi / th);
% x_l95 = (x_med - circ_confmean(th*x*pi/180,1-0.95)'* 180/pi / th);
% x_u68 = (x_med + circ_confmean(th*x*pi/180,1-0.68)'* 180/pi / th);
% x_l68 = (x_med - circ_confmean(th*x*pi/180,1-0.68)'* 180/pi / th);

x_u95 = (x_med + 2*circ_std(th*x*pi/180)'* 180/pi /th);
x_l95 = (x_med - 2*circ_std(th*x*pi/180)'* 180/pi /th);
x_u68 = (x_med + circ_std(th*x*pi/180)'* 180/pi /th);
x_l68 = (x_med - circ_std(th*x*pi/180)'* 180/pi /th);

% Make sure all angles are between -180 and 180
for ii = flip(1:5)
    x_med(x_med<(ii-1)*180) = x_med(x_med<(ii-1)*180) + ii*180;
    x_u95(x_u95<(ii-1)*180) = x_u95(x_u95<(ii-1)*180) + ii*180;
    x_l95(x_l95<(ii-1)*180) = x_l95(x_l95<(ii-1)*180) + ii*180;
    x_u68(x_u68<(ii-1)*180) = x_u68(x_u68<(ii-1)*180) + ii*180;
    x_l68(x_l68<(ii-1)*180) = x_l68(x_l68<(ii-1)*180) + ii*180;

    x_med(x_med>ii*180) = x_med(x_med>ii*180) - ii*180;
    x_u95(x_u95>ii*180) = x_u95(x_u95>ii*180) - ii*180;
    x_l95(x_l95>ii*180) = x_l95(x_l95>ii*180) - ii*180;
    x_u68(x_u68>ii*180) = x_u68(x_u68>ii*180) - ii*180;
    x_l68(x_l68>ii*180) = x_l68(x_l68>ii*180) - ii*180;
end


end

