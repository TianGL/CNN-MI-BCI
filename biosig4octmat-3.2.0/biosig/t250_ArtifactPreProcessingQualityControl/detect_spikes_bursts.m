function [HDR, s] = detect_spikes_bursts(fn, chan, varargin)
% DETECT_SPIKES_BURSTS detects spikes and bursts of spikes in
% neural recordings.
% Spikes are detected when voltages increase is larger than slopeThreshold
% (default 20 V/s) within a window of length winlen (default 0.0002 s). An interspike interval
% large than dT_Burst (default 75 ms) define the start of the next burst.
%
%
% HDR = detect_spikes_bursts(filename, chan)
% ... = detect_spikes_bursts(HDR, data)
% ... = detect_spikes_bursts(... ,'OptimumJK')
% ... = detect_spikes_bursts(... ,'-o',outputFilename)
% ... = detect_spikes_bursts(... ,'-e',eventFilename)
% ... = detect_spikes_bursts(... ,'-b',burstFilename)
% ... = detect_spikes_bursts(... ,'slopeThreshold',slopeThreshold)
% ... = detect_spikes_bursts(... ,'winlen',winlen)
% ... = detect_spikes_bursts(... ,'dT_Burst',dT_Burst)
% ... = detect_spikes_bursts(... ,'dT_Exclude',dT_Exclude)
% ... = detect_spikes_bursts(... ,[slopeThreshold [, winlen [, dT_Burst [,dT_Exclude ]]] ])
% ... = detect_spikes_bursts(... ,'OptimumJK')
% [HDR, data] = detect_spikes_bursts(...)
%
% Input:
% 	filename: name of source file
%	chan	list of channels that should be analyzed (default is 0: all channels)
%	HDR	header structure obtained by SOPEN, SLOAD, or meXSLOAD
%	data	signal data that should be analyzed
%	slopeThreshold	[default: 10 V/s] Spike is detected when
%		slope (over time winlen) exceeds this value
%	winlen	[default: .2e-3 s] windowlength in seconds for computing slope
%	dT_Burst	[default: 50e-3 s] am inter-spike-interval (ISI) exceeding this value,
%		marks the beginning of a new burst
%               'OPTIMUM_ISI' will use the function OPTIMUM_ISI_SPIKE_BURST_SEPARATION
%               for identifying dT_Burst.
%	dT_Exclude an interspike interval smaller than this value, indicates a
%		double detection, and the second detection is deleted.
%		in case of several consecutive ISI's smaller than this value,
%		all except the first spikes are deleted.
%	'OptimumJK' overrides previously defined settings, and uses the optimal setting based on
%		an JK's data set. It sets slopeThreshold = 2.5, dT=450e-6, dT_Exclude = 2e-3;
%		These settings were found based on 18 different recordings containing
%		23999 spikes, 3967 bursts, and over 9700 seconds of recorded data.
%	outputFilename
%		name of file for storing the resulting data with the
%		detected spikes and bursts in GDF format.
%	eventFilename
%		filename to store event inforamation in GDF format. this is similar to
%		the outputFile, except that the signal data is not included and is, therefore,
%		much smaller than the outputFile
%	burstFilename
%		filename for the "burst table", containing basic properties of each burst,
%		(it is an ASCII file in <tab>-delimited format)
%
%	Arguments can appear in any order and multiple times (except for filename, chan, HDR and data),
%	In case of conflicting definitions, the latest definition has highest precedence and is used.
%
%
% Output:
%     HDR	header structure as defined in biosig
%     HDR.EVENT includes the detected spikes and bursts.
%     HDR.BurstTable contains for each burst (each in a row) the following 5 numbers:
%	channel number, sweep number, OnsetTime within sweep [s],
%	number of spikes within burst, and average inter-spike interval (ISI) [ms]
%     data	signal data, one channel per column
%		between segments, NaN values for 0.1s are introduced
%
% see also:  DETECT_SPIKES_BURSTS, SPIKE2BURSTS, OPTIMUM_ISI_SPIKE_BURST_SEPARATION
%
% References:
%

%    Copyright (C) 2011,2014,2016 by Alois Schloegl <alois.schloegl@ist.ac.at>
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.

if nargin < 1,
	help detect_spikes_bursts;
end;

if nargin < 2,
	chan = 0;
end;

%%%%% default settings %%%%%
dT = .2e-3;		%%% smoothing window length [s]
dT_Burst = 50e-3;	%%% smaller ISI is a burst [s]
dT_Exclude = [];	%%% for identifying double detections, spikes with smaller ISI are excluded
slopeThreshold = 10; 	%% [V/s]
outFile = [];
evtFile = [];
burstFile = [];


%%%%% analyze input arguments %%%%%
k = 1;
K = 0;
while k <= length(varargin)
	if ischar(varargin{k})
		if (strcmp(varargin{k},'-o'))
			k = k + 1;
			outFile = varargin{k};
		elseif (strcmp(varargin{k},'-e'))
			k = k + 1;
			evtFile = varargin{k};
		elseif (strcmp(varargin{k},'-b'))
			k = k + 1;
			burstFile = varargin{k};
		elseif (strcmpi(varargin{k},'slopeThreshold'))
			k = k + 1;
			slopeThreshold = varargin{k};
		elseif (strcmpi(varargin{k},'winlen'))
			k = k + 1;
			dT = varargin{k};
		elseif (strcmpi(varargin{k},'dT_Burst'))
			k = k + 1;
			dT_Burst = varargin{k};
		elseif (strcmpi(varargin{k},'dT_Exclude'))
			k = k + 1;
			dT_Exclude = varargin{k};
		elseif (strcmpi(varargin{k},'OptimumJK'))
			slopeThreshold = 2.5; %
			dT = .450e-3;
			dT_Exclude = 2e-3;
		else
			warning(sprintf('unknown input argument <%s>- ignored',varargin{k}));
		end;
	elseif isnumeric(varargin{k})
		K = K + 1
		switch (K)
		case {1}
			slopeThreshold = varargin{k}; %% [V/s]
		case {2}
			dT = varargin{k};		%%% smoothing window length [s]
		case {3}
			dT_Burst = varargin{k};		%%% smaller ISI is a burst [s]
		case {4}
			dT_Exclude = varargin{k};	%%% smaller ISI is a burst [s]
		otherwise
			warning(sprintf('unknown input argument <%f> - ignored',varargin{k}));
		end;
	end;
	k = k+1;
end;


Fs = 20000; 	% assumed samplerate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ischar(fn) && exist(fn,'file') && any(size(chan)==1)
		winlen = ceil(Fs*.1);
		[s, HDR] = sload(fn, 0, 'NUMBER_OF_NAN_IN_BREAK', winlen);
		if Fs < HDR.SampleRate,
			winlen   = ceil(HDR.SampleRate * .1);
			[s, HDR] = sload(fn, 0, 'NUMBER_OF_NAN_IN_BREAK', winlen);
		end;
		if chan==0, chan = 1:HDR.NS; end;
	elseif isstruct(fn)
		HDR = fn;
		s = chan;
		HDR.NS = size(s,2);
		chan = 1:HDR.NS;
	else
		help(mfilename);
	end


	EVENT = HDR.EVENT;
	if ~isfield(EVENT,'DUR');
		EVENT.DUR = zeros(size(EVENT.POS));
	end;
	if ~isfield(EVENT,'CHN');
		EVENT.CHN = zeros(size(EVENT.TYP));
	end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Set Parameters for Spike Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	B  = [1; zeros(round(HDR.SampleRate*dT)-1, 1); -1];

	HDR.BurstTable = [];

	for ch = chan(:)';	% look for each channel
		% only voltage channels are considered
	if (bitand(HDR.PhysDimCode(ch), hex2dec('ffe0')) == 4256), %% physicalunits('V'),
		%%%%%%%	Spike Detection %%%%%%%%%%%%%%%%%%%
		[unit, scale] = physicalunits(HDR.PhysDimCode(ch));
		tmp = scale * filter(B, round(HDR.SampleRate*dT)/HDR.SampleRate, s(:,ch));
		OnsetSpike = find( diff (tmp > slopeThreshold) > 0);	%% spike onset time [samples]
		% --- remove double detections < 1 ms
		if ~isempty(dT_Exclude) && ~isempty(OnsetSpike),
			OnsetSpike = OnsetSpike([1; 1+find(diff(OnsetSpike) > Fs * dT_Exclude)]);
		end;


		EVENT.TYP = [EVENT.TYP; repmat(hex2dec('0201'), size(OnsetSpike))];
		EVENT.POS = [EVENT.POS; OnsetSpike];
		EVENT.DUR = [EVENT.DUR; repmat(0,  size(OnsetSpike))];
		EVENT.CHN = [EVENT.CHN; repmat(ch, size(OnsetSpike,1), 1) ];
		if isfield(EVENT,'TimeStamp')
			%%% TODO: these should be properly computed %%%
			EVENT.TimeStamp = [EVENT.TimeStamp; repmat(NaN, size(OnsetSpike,1), 1) ];
		end;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%	get peak time and max slope time
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if 1,
			POS = repmat(0, size(OnsetSpike).*[2,1]);
			TYP = repmat(NaN, size(POS));

			for k=1:length(OnsetSpike),
				% find peak time within interval of [-winlen ms, 0.001 ms]
				tix = OnsetSpike(k) + [0 : 0.001*HDR.SampleRate];
				tix = tix(tix<=size(s,1));
				%[peak,pix]=max(s(tix,ch)); pix=pix(1);
				stmp= s(tix,ch);
				pix = tix(find(stmp == max(stmp)));
				POS(2*k-1) = round(mean(pix));
				TYP(2*k-1) = hex2dec('204');

				% find max slope time
				tix = OnsetSpike(k) + [-round(HDR.SampleRate*dT) : 0.001*HDR.SampleRate];
				tix = tix(0<tix & tix<=size(tmp,1));
				stmp= tmp(tix);
				pix = tix(find(stmp == max(stmp)))-round(HDR.SampleRate*dT)/2;
				POS(2*k) = round(mean(pix));
				TYP(2*k) = hex2dec('203');
			end;
			EVENT.TYP = [EVENT.TYP; TYP];
			EVENT.POS = [EVENT.POS; POS];
			EVENT.DUR = [EVENT.DUR; repmat(0,  size(POS)) ];
			EVENT.CHN = [EVENT.CHN; repmat(ch, size(POS)) ];
			if isfield(EVENT,'TimeStamp')
				%%% TODO: these should be properly computed %%%
				EVENT.TimeStamp = [EVENT.TimeStamp; repmat(NaN, size(POS)) ];
			end;

		end
	end;
	end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Set Parameters for Burst Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	HDR.EVENT = EVENT;
	vararg = {'dT_Burst', dT_Burst, 'dT_Exclude', dT_Exclude};
	if ~isempty(evtFile)
		vararg{end+1}='-e';
		vararg{end+1}=evtFile;
	end;
	if ~isempty(burstFile)
		vararg{end+1}='-b';
		vararg{end+1}=burstFile;
	end;
	HDR = spikes2bursts(HDR,vararg{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	HDR.NRec = size(s,1);
	HDR.SPR = 1;
	HDR = rmfield(HDR,'AS');
	if ~isempty(outFile)
		%%% write data to output
		HDR.TYPE  = 'GDF';
		HDR.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		HDR.FILE = [];
		HDR.FileName  = outFile;
		HDR.FILE.Path = '';
		HDR.PhysMax   = max(s);
		HDR.PhysMin   = min(s);
		HDR.DigMax(:) =  2^15-1;
		HDR.DigMin(:) = -2^15;
		HDR.GDFTYP(:) = 3;
		HDR.FLAG.UCAL = 0;
		HDR.NS  = size(s,2);
		HDR.Dur = 1/HDR.SampleRate;
		HDR = sopen(HDR,'w');
		if (HDR.FILE.FID < 0)
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',HDR.FileName);
		else
			HDR = swrite(HDR,s);
			HDR = sclose(HDR);
		end;
	end;

