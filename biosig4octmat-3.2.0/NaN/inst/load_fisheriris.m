% LOAD_FISHERIRIS 
%        loads famous iris data set from Fisher, 1936 [1]. 
%
% References: 
% [1] Fisher,R.A. "The use of multiple measurements in taxonomic problems" 
%        Annual Eugenics, 7, Part II, 179-188 (1936); also in "Contributions to Mathematical Statistics" (John Wiley, NY, 1950).
% [2] Duda,R.O., & Hart,P.E. (1973) Pattern Classification and Scene Analysis. 
%        (Q327.D83) John Wiley & Sons. ISBN 0-471-22361-1. See page 218.

%	$Id$
%	Copyright (C) 2009,2010 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street - Fifth Floor, Boston, MA 02110-1301, USA.


if exist('OCTAVE_VERSION','builtin')
	if ~exist('iris.data','file')
		if strncmp(computer,'PCWIN',5)
			fprintf(1,'Download http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data and save in local directory %s\nPress any key to continue ...\n',pwd);
		else 
		        unix('wget http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data'); 
		end;         
	end;        
        tmp = fopen('iris.data'); species=fread(tmp,[1,inf],'uint8=>char'); fclose(tmp); 
        [meas,tmp,species]=str2double(species,',');        
        meas = meas(:,1:4);
        species = species(:,5);        
else 
        load fisheriris; 
end; 