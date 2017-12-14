function [h, pval, ci, stats] = ttest (x, m, alpha, tail, vartype, DIM)
% TTEST (paired) t-test
%     For a sample X from a normal distribution with unknown mean and
%     variance, perform a t-test of the null hypothesis `mean (X) == M'.
%     Under the null, the test statistic T follows a Student
%     distribution with `DF = length (X) - 1' degrees of freedom.
%
%     TTEST treads NaNs as "Missing values" and ignores these. 
%
% H = ttest(x,m)
%	tests Null-hypothesis that mean of x is m. 		
% H = ttest(x,y)
% 	size of x and size of y must match, it is tested whether the 
%	difference x-y is significantly different to m=0; 
% H = ttest(x,y,alpha)
% H = ttest(x,y,alpha,tail)
% H = ttest(x,y,alpha,tail,DIM)
% [H,PVAL] = ttest(...)
%
%     H=1 indicates a rejection of the Null-hypothesis at a significance 
%     level of alpha (default alpha = 0.05).	 
% 
%     With the optional argument string TAIL, the alternative of interest
%     can be selected.  If TAIL is '!=' or '<>' or 'both', the null is tested
%     against the two-sided Alternative `mean (X) ~= mean (Y)'.  If TAIL
%     is '>' or 'right', the one-sided Alternative `mean (X) > mean (Y)' is used.
%     Similarly for '<' or 'left', the one-sided Alternative `mean (X) < mean
%     (Y)' is used.  The default is the two-sided case.
% 
%     H returns whether the Null-Hypotheses must be rejected. 
%     The p-value of the test is returned in PVAL. 
% 
%     TTEST works on the first non-singleton dimension or on DIM. 
%
%     If no output argument is given, the p-value of the test is
%     displayed.
%

%%% not supported yet 
% [h,p,ci] = ttest(...)
% [h,p,ci,stats] = ttest(...)

%       $Id$
%       Copyright (C) 1995, 1996, 1997, 1998, 2000, 2002, 2005, 2006, 2007
%               Kurt Hornik
%       Copyright (C) 2010 by Alois Schloegl <alois.schloegl@gmail.com>	
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

  if ((nargin < 2) || (nargin > 5) || nargout > 4)
        print_usage ;
  end

  if (nargin == 2)
    alt  = '~=';
  end
  if (nargin < 3) || isempty(alpha)
    alpha = .05;
  end

  if (nargin < 4) || isempty(tail)
    tail = '~=';
  end
  if (~ ischar (tail))
    error ('ttest: tail must be a string');
  end
  if (nargin < 5) || isempty(vartype)
    vartype = 'equal';
  end
  if ~strcmp(vartype,'equal')
    error ('test: vartype not supported')
  end	
  if nargin<6,
       DIM = find(size(x)>1,1);
  end;
  if isempty(DIM), DIM=1; end;


  szx = size(x); 
  szm = size(m);	
  szx(DIM) = 1;	  
  szm(DIM) = 1;
  if size(m,DIM)==1
  	;
  elseif size(x,DIM) == size(m,DIM)
  	x = x-m;
  	m = zeros(szm);
  else
    error ('ttest: dimension of X and Y do not fit');
  end 	  
  
  [S, N] = sumskipnan(x, DIM);
  stats.df = N - 1;
  stats.sd = std (x);
  stats.tstat = sqrt (N) .* (S./N - m) ./ stats.sd;
  cdf = tcdf (stats.tstat, stats.df);

  if (strcmp (tail, '~=') || strcmp (tail, '!=') || strcmp (tail, '<>')) || strcmp(tail,'both'),
    pval = 2 * min (cdf, 1 - cdf);
  elseif strcmp (tail, '>') || strcmp(tail,'right'),
    pval = 1 - cdf;
  elseif strcmp (tail, '<') || strcmp(tail,'left'),
    pval = cdf;
  else
    error ('ttest: option %s not recognized', tail);
  end

  h = pval < alpha;	
  if (nargout == 0)
    fprintf(1,'  pval: %g\n', pval);
  end

