function [RES] = roc(d, c, varargin);
% ROC plots receiver operator curve and computes derived statistics.
%   In order to speed up the plotting, no more than 10000 data points are shown at most.
%   [if you need more, you need to change the source code]
%
% Usage:
% [...] = roc(d,c);
% [...] = roc(d1,d0);
% [...] = roc(...,s);
% [SEN, SPEC, TH, ACC, AUC,Yi,idx]=roc(...);
% RES = roc(...,'FPR',FPR);
%	RES.THRESHOLD.FPR returns the threshold value to obtain
%	the given FPR rate.
% RES = roc(...,'maxYI');
% RES = roc(...,'maxACC');
% RES = roc(...,'maxKAPPA');
%	RES.THRESHOLD.{maxYI,maxACC,maxKAPPA} return the threshold
%	value to obtain the maxium YoudenIndex (YI), Accuracy and Kappa, resp.
%
% INPUT:
% d	DATA
% c	CLASS, vector with 0 and 1
% d1	DATA of class 1
% d2	DATA of class 0
% s	line style (as used in plot)
%
% OUTPUT:
% SEN     sensitivity
% SPEC    specificity
% TH      Threshold
% ACC     accuracy
% AUC     area under ROC curve
% Yi	  max(SEN+SPEC-1), Youden index
% c	  TH(c) is the threshold that maximizes Yi
%
%   Remark: if the sample values in d are not unique, there is a certain
%      ambiguity in the results; the results may vary depending on
%      on the ordering of the samples. Usually, this is only a problem,
%      if the number of unique data value is much smaller than the total
%      number of samples.
%
% see also: AUC, PLOT, ROC

%	Copyright (c) 1997-2003,2005,2007,2010,2011,2016,2017 Alois Schloegl <alois.schloegl@gmail.com>
%	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.
%

global FLAG_DISPLAY;

MODE = all(size(d)==size(c)) && all(all((c==1) | (c==0) | isnan(c)));
d=d(:);
c=c(:);

if ~MODE
        d2=c;
        c=[ones(size(d));zeros(size(d2))];
        d=[d;d2];
        fprintf(2,'Warning ROC: XXX\n')
else
	ix = ~any(isnan([d,c]),2);
	c = c(ix);
	d = d(ix);
end;

% handle (ignore) NaN's
c = c(~isnan(d));
d = d(~isnan(d));

plot_args={'-'};
flag_plot_args = 1;
thFPR = NaN;
for k=1:length(varargin)
	arg = varargin{k};
	if strcmp(arg,'FPR')
		flag_plot_args = 0;
		thFPR = varargin{k+1};
	end;
	if flag_plot_args,
		plot_args{k} = arg;
	end
end;

[D,I] = sort(d);
x = c(I);

FNR = cumsum(x==1)/sum(x==1);
TPR = 1-FNR;

TNR = cumsum(x==0)/sum(x==0);
FPR = 1-TNR;

FN = cumsum(x==1);
TP = sum(x==1)-FN;

TN = cumsum(x==0);
FP = sum(x==0)-TN;

RES.PPV = TP./(TP+FP);
RES.NPV = TN./(TN+FN);

SEN = TP./(TP+FN);
SPEC= TN./(TN+FP);
ACC = (TP+TN)./(TP+TN+FP+FN);

% SEN = [FN TP TN FP SEN SPEC ACC D];

%%% compute Cohen's kappa coefficient
N = size(d,1);

%H = [TP,FN;FP,TN];
p_i = [TP+FP,FN+TN];%sum(H,1);
pi_ = [TP+FN,FP+TN];%sum(H,2)';
pe  = sum(p_i.*pi_,2)/(N*N);  % estimate of change agreement
kap = (ACC - pe) ./ (1 - pe);
mcc = (TP .* TN - FN .* FP) ./ sqrt(prod( [p_i, pi_], 2));

% mutual information
pxi = pi_/N;                       % p(x_i)
pyj = p_i/N;                       % p(y_j)
log2pji = ([TP,FN,FP,TN]/N).*log2([TP,FN,FP,TN]./[p_i,p_i]);
RES.MI = -sumskipnan(pyj.*log2(pyj),2) + sumskipnan(log2pji,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  display only 10000 points at most.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FLAG_DISPLAY,
	len = length(FPR);
	delta = max(1,floor(len/5000));
	ix = [1:delta:len-1,len];
	plot(FPR(ix)*100,TPR(ix)*100, plot_args{:});

	%ylabel('Sensitivity (true positive ratio) [%]');
	%xlabel('1-Specificity (false positive ratio) [%]');
end;

% area under the ROC curve
RES.AUC = -diff(FPR)' * (TPR(1:end-1)+TPR(2:end))/2;

% Youden index
YI = SEN + SPEC - 1;

RES.YI    = YI;
RES.ACC   = ACC;
RES.KAPPA = kap;
RES.MCC   = mcc;
RES.TH    = D;
RES.F1    = 2*TP./(2*TP+FP+FN);

RES.SEN = SEN;
RES.SPEC = SPEC;
RES.TP = TP;
RES.FP = FP;
RES.FN = FN;
RES.TN = TN;
RES.TPR = TPR;
RES.FPR = FPR;
RES.FNR = FNR;
RES.TNR = TNR;
RES.LRP = TPR./FPR;
RES.LRN = FNR./TNR;

% find optimal threshold
[tmp,ix] = max(SEN+SPEC-1);
RES.THRESHOLD.maxYI = D(ix);
RES.H_yi = [TP(ix),FN(ix);FP(ix),TN(ix)];

[RES.maxKAPPA,ix] = max(kap);
RES.THRESHOLD.maxKAPPA = D(ix);
RES.H_kappa = [TP(ix),FN(ix);FP(ix),TN(ix)];

[RES.maxMCC,ix] = max(mcc);
RES.THRESHOLD.maxMCC = D(ix);
RES.H_mcc = [TP(ix),FN(ix);FP(ix),TN(ix)];

[RES.maxMI,ix] = max(RES.MI);
RES.THRESHOLD.maxMI = D(ix);
RES.H_mi = [TP(ix),FN(ix);FP(ix),TN(ix)];

[tmp,ix] = max(ACC);
RES.THRESHOLD.maxACC = D(ix);
RES.H_acc = [TP(ix),FN(ix);FP(ix),TN(ix)];

[tmp,ix] = max(RES.F1);
RES.THRESHOLD.maxF1 = D(ix);
RES.H_f1 = [TP(ix),FN(ix);FP(ix),TN(ix)];

RES.THRESHOLD.FPR = NaN;
if ~isnan(thFPR)
	ix = max(1,min(N,round((1-thFPR)*N)));
	RES.THRESHOLD.FPR = D(ix);
	RES.H_fpr = [TP(ix),FN(ix);FP(ix),TN(ix)];
end;
