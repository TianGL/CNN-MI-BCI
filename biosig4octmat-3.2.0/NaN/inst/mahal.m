function [d] = mahal(X,Y)
% MAHAL return the Mahalanobis' D-square distance between the 
% multivariate samples x and y, which must have the same number 
% of components (columns), but may have a different number of observations (rows). 
% 
%  d = mahal(X,Y)
%
%   d(k) = (X(k,:)-MU)*inv(SIGMA)*(X(k,:)-MU)'
%
%   where MU and SIGMA are the mean and the covariance matrix of Y 
%
%
% see also: TRAIN_SC, TEST_SC, COVM
%
% References: 

%	$Id: train_sc.m 2140 2009-07-02 12:03:55Z schloegl $
%	Copyright (C) 2009 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
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

sx = size(X);
sy = size(Y);

if sx(2)~=sy(2),
	error('number of columns of X and Y do not fit');
end;	

IR = inv(covm(Y,'E'));
IR(1,1) = IR(1,1)-1;

D = [ones(size(X,1),1),X];
d = sum((D*IR).*D,2);

