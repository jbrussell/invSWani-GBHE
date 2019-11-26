function [ x_med, x_u95, x_l95, x_u68, x_l68 ] = get_95_68_prctile( x )
% [ x_med, x_u95, x_l95, x_u68, x_l68 ] = get_95_68_prctile( x )
%
% Get the upper and low 95th and 68th percentiles and median;
%
x_u95 = prctile(x,97.5,1);
x_l95 = prctile(x,2.5,1);
x_u68 = prctile(x,16,1);
x_l68 = prctile(x,84,1);
x_med = prctile(x,50,1);

end

