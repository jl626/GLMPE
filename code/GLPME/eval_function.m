function curr_norm = eval_function(fg,Fp,x)
%EVAL_FUNCTION internal evaluation function
%
%   Url: https://lts2.epfl.ch/unlocbox/doc/solver/misc/eval_function.php

% Copyright (C) 2012-2016 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.7.3
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


curr_norm = fg.eval(x); 


for ii = 1:length(Fp)
        curr_norm = curr_norm + Fp{ii}.eval(x); 
end

end
