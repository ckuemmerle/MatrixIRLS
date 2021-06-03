function in = inside(x, a, b)
% call: in = inside(x, a, b)
% input: 3 scalars x, a, b
% output: 1 (true) if x lies between a and b
%         0 (false) if x does not lie between a and b
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

% this will work even if x is infinity or NaN
% but first, make sure x is real
in = 0;
if ~isreal(x)   
   return
end
if a <= b
   if x >= a & x <= b
      in = 1;
   end
else
   if x >= b & x <= a
      in = 1;
   end
end
