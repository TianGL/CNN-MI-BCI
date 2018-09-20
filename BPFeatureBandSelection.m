function [muBand, betaBand] = BPFeatureBandSelection(data, labels, fs)
%****************************************************************
% subject optimal frequency selection based on BP FDA-F-Score
%
% [muBand, betaBand] = FeatureBandSelection[C3, C4, Cz]
% selection subject optimal frequency band from electrodes of C3, CZ and Cz
% data: eeg data of C3, C4, Cz electrode,
% fs: sample frequency
% [muBand, betaBand]:optimal mu band and beta band
% Author: Geliang Tian (tglasd@163.com)
% Copyright 2018
%*****************************************************************

% initialization parameter
muStep = 1; betaStep = 1;
muRange = [4,14]; betaRange = [16,40];
muSize = [4,5,6]; betaSize = [10,12];


size = length(data);

%%
% mu
nF = 0; % flag of frequency band
for kk = 1:length(muSize)
    for n_mu = 0 : floor( (muRange(2)-muRange(1)-muSize(kk)) / muStep )
        
        % 6th (5th may better) Butterworth filter
        mu0 = [muRange(1)+n_mu*muStep, muRange(1) + muSize(kk) + n_mu*muStep]
        
        % filter design
        muFilter = design( fdesign.bandpass('N,F3dB1,F3dB2', 6, mu0(1), mu0(2), fs), 'butter');
        n1 = 0;
        n2 = 0;
        
        % caculate feature value
        for i = 1:size
            muData = filter(muFilter, data{i} );
            if(labels(i) == 1)
                n1 = n1 + 1;
                BP1.mu(n1, :) = log10( var( muData ));
            else
                n2 = n2 + 1;
                BP2.mu(n2, :) = log10( var( muData ));
            end
        end
        BP1.mean = mean(BP1.mu);
        BP2.mean = mean(BP2.mu);
        
        S1.mu = var( BP1.mu );
        S2.mu = var( BP2.mu);
        % caculate F-score
        nF = nF + 1;
        Res.FscoreMu(nF) = sum((BP1.mean - BP2.mean).^2) / sum(S1.mu + S2.mu);
        Res.mu(nF,:) = mu0;
    end
end
%%
%beta
nF = 0; % flag of frequency band
for kk = 1:length(betaSize)
    for n_beta = 0 : floor( (betaRange(2)-betaRange(1)-betaSize(kk)) / betaStep )
        
        % 6th Butterworth filter
        beta0 = [betaRange(1)+n_beta*betaStep, betaRange(1) + betaSize(kk) + n_beta*betaStep]
         % filter design
        betaFilter = design( fdesign.bandpass('N,F3dB1,F3dB2', 6, beta0(1), beta0(2), fs), 'butter');
        
        n1 = 0;
        n2 = 0;
        
        % caculate feature value
        for i = 1:size
            betaData= filter(betaFilter, data{i} );
            if(labels(i) == 1)
                n1 = n1 + 1;
                BP1.beta(n1, :) = log10( var( betaData ));
            else
                n2 = n2 + 1;
                BP2.beta(n2, :) = log10( var( betaData ) );
            end
        end
        
        BP1.mean = mean(BP1.beta);
        BP2.mean = mean(BP2.beta);
        
        S1.beta = var( BP1.beta );
        S2.beta = var( BP2.beta);
        % caculate F-score
        nF = nF + 1;
        Res.FscoreBeta(nF) = sum((BP1.mean - BP2.mean).^2) / sum(S1.beta + S2.beta);
        Res.beta(nF,:) = beta0;
    end
end
num = find(Res.FscoreMu == max(Res.FscoreMu), 1, 'first');
muBand = Res.mu(num,:);
num = find(Res.FscoreBeta == max(Res.FscoreBeta), 1, 'first');
betaBand = Res.beta(num,:);
end