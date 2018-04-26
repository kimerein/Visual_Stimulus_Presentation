function map = jetm(n)
%JETM              Muted jet colormap.
%   JETM(N), a variant of JET(N), is an N-by-3 matrix containing a muted
%   version of the JET colormap.  JETM with no arguments is the same
%   length as the current colormap.

if (nargin < 1), n = size(colormap,1);  end;

% This colormap was constructed manually by altering the JET colormap,
%  which is visually a little too bright and suffers from light blue and
%  yellow bands.  To obtain a muted version,
%    (1) map a jet colormap to HSV color space
%    (2) linearly rescale the saturation values from the
%               range [0.5 1.0] => [0.25 0.65] (softens the colors)
%    (3) linearly rescale the color value numbers from the
%               range [0.5 1.0] => [0.8 1.0]  (keeps the colors bright)
%    (4) use the colormapeditor to better space out the light blue/yellow
%               color points
%    (6) remap from HSV to RGB color space
%    (5) symmetrize the map by enforcing the red/blue values to be
%               inverted copies of one another and the green values
%               to be symmetric (evens out the colormapeditor step)
% The result (sampled at 64 color points) is:

map = ...
[...
0.2800 0.2811 0.8062; ...
0.2912 0.2952 0.8547; ...
0.3024 0.3092 0.9031; ...
0.3136 0.3233 0.9516; ...
0.3248 0.3374 1.0000; ...
0.3360 0.3688 0.9946; ...
0.3472 0.4002 0.9893; ...
0.3508 0.4493 0.9839; ...
0.3519 0.5042 0.9785; ...
0.3530 0.5591 0.9731; ...
0.3541 0.6140 0.9678; ...
0.3552 0.6690 0.9624; ...
0.3563 0.7239 0.9570; ...
0.3519 0.7713 0.9516; ...
0.3419 0.7960 0.9494; ...
0.3320 0.8207 0.9472; ...
0.3220 0.8454 0.9450; ...
0.3121 0.8700 0.9428; ...
0.3021 0.8947 0.9406; ...
0.2922 0.9194 0.9384; ...
0.2837 0.9365 0.9290; ...
0.2797 0.9419 0.9126; ...
0.2757 0.9472 0.8962; ...
0.2717 0.9526 0.8797; ...
0.2677 0.9578 0.8589; ...
0.3005 0.9566 0.8251; ...
0.3456 0.9533 0.7912; ...
0.3907 0.9501 0.7574; ...
0.4358 0.9468 0.7235; ...
0.4809 0.9435 0.6896; ...
0.5260 0.9402 0.6558; ...
0.5627 0.9394 0.6219; ...
0.5965 0.9394 0.5880; ...
0.6304 0.9394 0.5542; ...
0.6642 0.9411 0.5147; ...
0.6981 0.9443 0.4696; ...
0.7320 0.9476 0.4245; ...
0.7658 0.9509 0.3794; ...
0.7997 0.9542 0.3343; ...
0.8336 0.9574 0.2892; ...
0.8674 0.9566 0.2687; ...
0.8838 0.9512 0.2727; ...
0.9003 0.9459 0.2767; ...
0.9167 0.9406 0.2807; ...
0.9331 0.9352 0.2847; ...
0.9389 0.9133 0.2947; ...
0.9411 0.8886 0.3046; ...
0.9433 0.8639 0.3146; ...
0.9456 0.8392 0.3245; ...
0.9478 0.8145 0.3345; ...
0.9500 0.7898 0.3444; ...
0.9530 0.7613 0.3544; ...
0.9584 0.7101 0.3560; ...
0.9637 0.6552 0.3549; ...
0.9691 0.6003 0.3538; ...
0.9745 0.5454 0.3527; ...
0.9798 0.4905 0.3516; ...
0.9852 0.4356 0.3505; ...
0.9906 0.3924 0.3444; ...
0.9960 0.3610 0.3332; ...
0.9879 0.3339 0.3220; ...
0.9395 0.3198 0.3108; ...
0.8910 0.3057 0.2996; ...
0.8426 0.2917 0.2884; ...
];

map = interp1(linspace(0,1,size(map,1)), map, linspace(0,1,n), 'linear');