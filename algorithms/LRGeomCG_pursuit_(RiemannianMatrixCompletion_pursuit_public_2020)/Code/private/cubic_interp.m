function xmin = cubic_interp(x1, x2, f1, f2, g1, g2) 
% call: xmin = cubic_interp(x1, x2, f1, f2, g1, g2) 
% find minimizer of the cubic polynomial Hermite-interpolating a
% function of one variable, at the two points x1 and x2, using the
% function (f1 and f2) and derivative (g1 and g2) data.
% could easily return inf or NaN, so this must be checked after call

% The formula is given by Nocedal and Wright, but the one given in
% the first edition has the wrong sign when x1 > x2. 
% Found the correct formula in a 1974 NPL report by Gill and Murray,
% and it is corrected in the second edition of Nocedal and Wright.
% Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
% with a subject header containing the string "nlcg".
% NLCG Version 1.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  NLCG 1.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta = g1 + g2 - 3*(f1 - f2)/(x1 - x2);    % called d1 in Nocedal and Wright
gamma = sign(x2-x1)*sqrt(eta^2 - g1*g2);  % called d2 in Nocedal and Wright
xmin = x2 - (x2 - x1)*(g2 + gamma - eta)/(g2 - g1 + 2*gamma);
   


