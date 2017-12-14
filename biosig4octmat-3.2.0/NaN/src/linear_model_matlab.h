/*

$Id$
Copyright (c) 2007-2009 The LIBLINEAR Project.
Copyright (c) 2010 Alois Schloegl <alois.schloegl@gmail.com>
This function is part of the NaN-toolbox
http://pub.ist.ac.at/~schloegl/matlab/NaN/

This code was extracted from liblinear-1.51 in Jan 2010 and 
modified for the use with Octave 

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

*/


#ifdef __cplusplus
extern "C" {
#endif

const char *model_to_matlab_structure(mxArray *plhs[], struct model *model_);
const char *matlab_matrix_to_model(struct model *model_, const mxArray *matlab_struct);

#ifdef __cplusplus
}
#endif

