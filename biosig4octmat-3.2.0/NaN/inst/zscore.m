function [i,v,m] = zscore(i,DIM)
% ZSCORE removes the mean and normalizes data 
% to a variance of 1. Can be used for pre-whitening of data, too. 
%
% [z,r,m] = zscore(x,DIM)
%   z   z-score of x along dimension DIM
%   r   is the inverse of the standard deviation
%   m   is the mean of x
%
% The data x can be reconstucted with 
%     x = z*diag(1./r) + repmat(m,size(z)./size(m))  
%     z = x*diag(r) - repmat(m.*v,size(z)./size(m))  
%
% DIM	dimension
%	1: STATS of columns
%	2: STATS of rows
%	default or []: first DIMENSION, with more than 1 element
%
% see also: SUMSKIPNAN, MEAN, STD, DETREND
%
% REFERENCE(S):
% [1] http://mathworld.wolfram.com/z-Score.html

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; If not, see <http://www.gnu.org/licenses/>.


%	$Id: zscore.m 10929 2012-08-29 19:31:12Z nir-krakauer $
%	Copyright (C) 2000-2003,2009 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/


if any(size(i)==0); return; end;

if nargin<2
        DIM=[]; 
end
if nargin<3
        W = []; 
end
if isempty(DIM), 
        DIM=min(find(size(i)>1));
        if isempty(DIM), DIM=1; end;
end;


% pre-whitening
m = mean(i,DIM,W);
i = i-repmat(m,size(i)./size(m));  % remove mean
v = 1./sqrt(mean(i.^2,DIM,W));
i = i.*repmat(v,size(i)./size(v)); % scale to var=1

        
