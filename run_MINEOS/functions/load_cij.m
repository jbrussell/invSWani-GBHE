function [ cij ] = load_cij( )
% UNITS = Mbar*100 [ends up with velocity in units of km/s]
%         For units of m/s use Pa = Mbar*100*1e9
%
% 1 : E&D PYROLITE 87
% 2 : Pesilnick and Nicolas, 78 (Mesozoic)
% 3 : Estey & Douglas 86
% 4 : Ismail & Mainprice, 98 (fast-spreading ridge) [No rotation]
% 5 : Ismail & Mainprice, 98 (fast-spreading ridge) [Jim rotation]
% 6 : Estey & Douglas 86 (Hexagonal)
% 7 : A-type (Jung et al 2006)
% 8 : B-type (Jung et al 2006)
% 9 : C-type (Jung et al 2006)
% 10: E-type (Jung et al 2006)
% 11: Ismail & Mainprice, 98 (fast-spreading ridge) [Correct rotation]
% 12: NoMelt 35 km [ignoring directional terms]
% 13: Nomelt Moho (~11 km) [ignoring directional terms]
% 14: NoMelt 35 km [full tensor, rotated]
% 15: NoMelt Moho (~11 km) [full tensor, rotated]
% 16: E-type (Katayama et al 2004)
% 17: Kaapvaal (Ben-Ismail et al 2001)
% 18: B-type Imono Dunites (Tasaka et al 2008)

%% PYROLITE (horizontal)
% E&D PYROLITE 87
cij(1).ref = 'E&D PYROLITE 87';
cij(1).c = [269. 55. 63. 0. 0. 0.;
      55. 200. 67. 0. 0. 0.;
     63. 67. 197.5 0. 0. 0.;
         0. 0. 0. 70. 0. 0.;
         0. 0. 0. 0. 78. 0.;
         0. 0. 0. 0. 0. 80.];

%% [Pesilnick and Nicolas, 1978]
% % %Pesilnick and Nicolas, 1978] (Mesozoic) (THIS IS THE ONE THAT WORKS)
% cij(2).ref = 'P&N 78 (Mesozoic)';
% cij(2).c = [236.5 72.5 72.3  0.   0.   0.;
%     72.5  220.8 71.9  0.   0.   0.;
%    72.3   71.9 220.2 0.   0.   0.;
%     0.   0.   0.    74.9   0.   0.;
%     0.   0.   0.     0    79.2  0.;
%     0.   0.   0.     0     0.  78.8 ];   
% %[ c ] = MS_rotEuler( c, FSD, 0, 0, 'sense', 'active' ); % rotate fast direction to FSD
% % cij(2).c = MS_rot3(cij(2).c,90,0,0);

% % %Pesilnick and Nicolas, 1978] (Mesozoic) FULL TENSOR
cij(2).ref = 'P&N 78 (Mesozoic)';
cij(2).c = [236.5   72.5   72.3   -0.08  -1.96      0;
             72.5  220.8   71.9    1.69   1.64  -0.40;
             72.3   71.9  220.2    1.82  -0.24  -0.41;
            -0.08   1.69   1.82    74.9  -0.85  -1.03;
            -1.96   1.64  -0.24   -0.85   79.2   1.28;
                0  -0.40  -0.41   -1.03   1.28   78.8 ];   
%[ c ] = MS_rotEuler( c, FSD, 0, 0, 'sense', 'active' ); % rotate fast direction to FSD
% cij(2).c = MS_rot3(cij(2).c,90,0,0);
%% % %  TEST FOR CELIA
% Estey & Douglas 86; 
% (Montagner & Anderson, 89)
cij(3).ref = 'E&D 86';
cij(3).c = [323.7  71.6   66.4  0.      0.     0.;
           71.6   235.1  75.6  0.      0.     0.;
           66.4   75.6  197.6  0.      0.     0.;
             0.     0.     0.   64.62    0.     0.; 
             0.     0.     0.    0.     79.04   0.;
             0.     0.     0.    0.      0.   78.05];
 
% % Estey & Douglas 86; (MADE HEXAGONAL)
% % (Montagner & Anderson, 89)
% cij(6).ref = 'E&D 86 (hexagonal)';
% cij(6).c = [235.1  71.6   66.4  0.      0.     0.;
%            71.6   235.1  66.4  0.      0.     0.;
%            66.4   66.4  197.6  0.      0.     0.;
%              0.     0.     0.   64.62    0.     0.;
%              0.     0.     0.    0.     64.62   0.;
%              0.     0.     0.    0.      0.   (235.1-71.6)/2];

% Estey & Douglas 86; (MADE HEXAGONAL)
% (Montagner & Anderson, 89)
cij(6).ref = 'E&D 86 (hexagonal)';
cij(6).c = [323.7  71.6   71.6  0.      0.     0.;
           71.6   235.1  66.4  0.      0.     0.;
           71.6   66.4  235.1  0.      0.     0.;
             0.     0.     0.   (235.1-66.4)/2    0.     0.;
             0.     0.     0.    0.     64.62   0.;
             0.     0.     0.    0.      0.   64.62];
         
% Hexagonal in our coordinate system (where G fast is in 1 direction)        
%           d c c . . .
%           c a b . . .
%           c b a . . .
%           . . . x . .          x = (a-b)/2
%           . . . . e .
%           . . . . . e

% In the paper...
%           a b c . . .
%           b a c . . .
%           c c d . . .
%           . . . e . .          x = (a-b)/2
%           . . . . e .
%           . . . . . x
      
 
%% Ismail and Mainprice, Tectonophys, 1998. (C and G have same direction instead of 45 degrees apart???)
% % X -> parallel to lineation [1 0 0] (maximum p-wave velocity)
% % Z -> perpendicular to foliation [0 1 0] (min p-wave velocity)
% % Y -> structural direction (foliation)? [0 0 1]
% % original -- fast-spreading ridge (JOSH)
% cij(4).ref = 'I&M 98 (fast-spreading unrotated)';
% cij(4).c =[1.9591  0.7112  0.7221  0.     0.     0.;
%    0.7112  2.3936  0.7187  0.     0.     0.;
%    0.7221  0.7187  2.0372  0.     0.     0.;
%     0.     0.      0.      0.7105 0.     0.;
%     0.     0.      0.      0      0.6249 0.; 
%     0.     0.      0.      0      0.     0.7024]*100.;
% 
% % double-rotated -- fast spreading ridge -- FIXED JIM VERSION
cij(5).ref = 'I&M 98 (fast-spreading) [Jim]';
cij(5).c =[2.3936 0.7187  0.7221  0.     0.     0.;
  0.7187  2.0372  0.7112  0.     0.     0.;
  0.7221  0.7112  1.9591  0.     0.     0.;
  0.     0.      0.      0.6249 0.     0.;
  0.     0.      0.      0      0.7105 0.; 
  0.     0.      0.      0      0.     0.7024]*100.;
% cij(5).c = MS_rot3(cij(5).c,90,0,0);
% 
% cij(11).ref = 'I&M 98 (fast-spreading)';
% rot3ax =  MS_rot3(cij(4).c,0,0,90);
% cij(11).c = MS_rot3(rot3ax,90,0,0);
% % cij(11).c = MS_rot3(cij(4).c,0,0,90);

%% Ismail and Mainprice, Tectonophys, 1998. FULL TENSOR
% X -> parallel to lineation [1 0 0] (maximum p-wave velocity)
% Z -> perpendicular to foliation [0 1 0] (min p-wave velocity)
% Y -> structural direction (foliation)? [0 0 1]
% original -- fast-spreading ridge (JOSH)
cij(4).ref = 'I&M 98 (fast-spreading unrotated)';
cij(4).c =[195.91   71.12   72.21  -0.07  -0.24  -0.14;
           71.12   239.36   71.87   0.03  -0.36  -0.16;
           72.21    71.87  203.72   0.40  -0.39  -0.02;
           -0.07     0.03    0.40  71.05  -0.15  -0.60;
           -0.24    -0.36   -0.39  -0.15  62.49   0.01; 
           -0.14    -0.16   -0.02  -0.60   0.01  70.24];

cij(11).ref = 'I&M 98 (fast-spreading)';

% fast a-axis in 1 direction, slow b-axis in 3 direction (vertical)
cij(11).c =  MS_rot3(cij(4).c,0,0,90);
cij(11).c = MS_rot3(cij(11).c,90,0,0);

% % fast a-axis in 1 direction, slow b-axis in 2 direction
% cij(11).c = MS_rot3(cij(4).c,0,0,90);

%% Cij Table 2 from Jung et al 2006 (5 GPa, 1573 K)

% A-type
cij(7).ref = 'A-type (Jung et al 2006)';
cij(7).c = [236.3   84.5    81.5     0.4     3.4     0.3;
             84.5  218.5    82.9    -1.8     1.2     0.3;
             81.5   82.9   208.0    -1.3     6.1     0.2;
              0.4   -1.8    -1.3    64.9    -0.1    -1.9;
              3.4    1.2     6.1    -0.1    68.7     0.3;
              0.3    0.3     0.2    -1.9     0.3    66.6];
% cij(7).c = MS_rot3(cij(7).c,-90,0,0);

% cij(7).c = MS_axes(cij(7).c);
% cij(7).c = MS_rot3(cij(7).c,0,90,0);
% cij(7).c = MS_rot3(cij(7).c,90,0,0);

% B-type
cij(8).ref = 'B-type (Jung et al 2006)';
cij(8).c = [221.3   84.3    83.3    -0.3     0.8     1.6;
             84.3  223.5    81.7    -1.4     1.0     1.6;
             83.3   81.7   215.5    -1.3     1.1    -0.4;
             -0.3   -1.4    -1.3    68.9     0.7    -0.4;
              0.8    1.0     1.1     0.7    67.4    -0.4;
              1.6    1.6    -0.4    -0.4    -0.4    69.4];
        
% C-type
cij(9).ref = 'C-type (Jung et al 2006)';
cij(9).c = [223.2   83.6    83.3     0.3    -3.7     0.3;
             83.6  209.8    81.9     0.8     1.5     0.3;
             83.3   81.9   228.5     0.4    -5.9     0.2;
              0.3    0.8     0.4    67.9     0.2    -1.9;
             -3.7    1.5    -5.9     0.2    71.1     0.3;
              0.3    0.3     0.2    -1.9     0.3    66.6];
        
% E-type
cij(10).ref = 'E-type (Jung et al 2006)';
cij(10).c = [236.8   82.3    84.1    -0.6     0.4     0.1;
              82.3  207.7    82.7    -2.6    -0.3    -1.0;
              84.1   82.7   217.4    -2.1    -1.9    -0.7;
              -0.6   -2.6    -2.1    65.0     0.1     0.4;
               0.4   -0.3    -1.9     0.1    71.1    -1.4;
               0.1   -1.0    -0.7     0.4    -1.4    68.5];
% cij(10).c = MS_rot3(cij(10).c,90,0,0);
           
%% NoMelt
% 35 km Depth
cij(12).ref = 'NoMelt 35km old';
cij(12).c = [270.3277  101.1694   99.6442         0         0         0;
             101.1694  232.8870   98.3383         0         0         0;
              99.6442   98.3383  240.5028         0         0         0;
                    0         0         0   66.4387         0         0;
                    0         0         0         0   74.9134         0;
                    0         0         0         0         0   71.5947];
                
% cij(12).c = MS_rot3(cij(12).c,90,0,0);

% Moho
cij(13).ref = 'NoMelt Moho old';
cij(13).c = [241.3361   76.2117   75.5783         0         0         0;
              76.2117  216.2073   74.8440         0         0         0;
              75.5783   74.8440  218.5703         0         0         0;
                    0         0         0   68.5739         0         0;
                    0         0         0         0   74.9277         0;
                    0         0         0         0         0   72.7660];
                
%% NoMelt Full Calculations
% 35 km Depth
cij(14).ref = 'NoMelt 35km';
cij(14).c = [233.5486  101.9677   98.3893         0         0    2.1218;
             101.9677  268.0697   99.5932         0         0    5.1256;
              98.3893   99.5932  240.5028         0         0    0.2528;
                    0         0         0   74.5830    1.6405         0;
                    0         0         0    1.6405   66.7691         0;
               2.1218    5.1256    0.2528         0         0   72.3929];
% Rotate by FSD (78) to get a-axis aligned with 1 direction
cij(14).c = MS_rot3(cij(14).c,0,0,78);
% % Rotate to be like old version with b-axis horizontal (E-type)
% cij(14).c = MS_axes(cij(14).c);
% cij(14).c = MS_rot3(cij(14).c,0,90,0);
% cij(14).c = MS_rot3(cij(14).c,90,0,0);

% Moho
cij(15).ref = 'NoMelt Moho';
cij(15).c = [217.8307   76.7795   74.9078         0         0    2.2513;
              76.7795  238.5773   75.5145         0         0    4.8380;
              74.9078   75.5145  218.5703         0         0    0.2069;
                    0         0         0   74.3737    1.7924         0;
                    0         0         0    1.7924   69.1278         0;
               2.2513    4.8380    0.2069         0         0   73.3337];
% Rotate by FSD (78) to get a-axis aligned with 1 direction
cij(15).c = MS_rot3(cij(15).c,0,0,78);
% % Rotate to be like old version with b-axis horizontal (E-type)
% cij(15).c = MS_axes(cij(15).c);
% cij(15).c = MS_rot3(cij(15).c,0,90,0);
% cij(15).c = MS_rot3(cij(15).c,90,0,0);

%% Katayama et al. 2004 (E-type experiments)
% E-type
cij(16).ref = 'E-type (Katayama et al 2004)';
cij(16).c = [258.1   79.0    87.3    -0.8   -13.8    -5.5;
              79.0  194.4    82.9    -1.7     2.1    -3.1;
              87.3   82.9   216.6    -1.9    -5.4    -1.3;
              -0.8   -1.7    -1.9    60.9    -1.8    -2.9;
             -13.8    2.1    -5.4    -1.8    76.2     0.4;
              -5.5   -3.1    -1.3    -2.9     0.4    67.8];
          
%% Ben Ismail, Barruol and Mainprice (GRL2001) BBM
% average of 48 Kaapvaal aggregates
%
cij(17).ref = 'Kaapvaal (BI et al 2001)';
cij(17).c =[2.2270  0.7855  0.7998  0.0009 -0.0010 -0.0003;
  0.7855  2.3639  0.7979 -0.0068 -0.0005  0.0031;
  0.7998  0.7979  2.2718 -0.0058  0.0009  0.0015;
  0.0009 -0.0068 -0.0058  0.7361  0.0018 -0.0008;
 -0.0010 -0.0005  0.0009  0.0018  0.7115 -0.0022;
 -0.0003  0.0031  0.0015 -0.0008 -0.0022  0.7494].*100;
cij(17).c =  MS_rot3(cij(17).c,0,0,90);
cij(17).c = MS_rot3(cij(17).c,90,0,0);

%% B-type Imono Dunites from fore-arc mantle wedge (Tasaka et al 2008)
%
% WARNING :: There's a typo in the Cijs provided! c12 is same as c11
% and pole plots look nothing like in the paper. I've changed c12 such that
% the pole plots look more like those from the paper.
%
% Also, it is also rotated such that [100] is along x3 and [010] is along x1.
% Instead, B-type should have [100] along x2 and [010] along x3, so must
% double rotate...

% cij(18).ref = 'B-type Imono Dunites (Tasaka 2008) [TYPO]';
% cij(18).c = [213.23  213.23    78.71  -0.0017    -0.09     -1.2;
%              213.23  226.38    81.11  -0.0042    -0.01    -1.39;
%               78.71   81.11   232.95     -0.2     0.21    -0.26;
%             -0.0017 -0.0042     -0.2    73.65    -0.46     0.04;
%               -0.09   -0.01     0.21    -0.46    71.55    -0.13;
%                -1.2   -1.39    -0.26     0.04    -0.13    69.15];
           
cij(18).ref = 'B-type Imono Dunites (Tasaka 2008)';
cij(18).c = [213.23   82.00    78.71  -0.0017    -0.09     -1.2;
              82.00  226.38    81.11  -0.0042    -0.01    -1.39;
              78.71   81.11   232.95     -0.2     0.21    -0.26;
            -0.0017 -0.0042     -0.2    73.65    -0.46     0.04;
              -0.09   -0.01     0.21    -0.46    71.55    -0.13;
               -1.2   -1.39    -0.26     0.04    -0.13    69.15];
cij(18).c =  MS_rot3(cij(18).c,90,0,0);
cij(18).c = MS_rot3(cij(18).c,0,90,0);

         
end

