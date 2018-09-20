function [muBand, betaBand] = ARFeatureBandSelection(data, labels, fs)
%****************************************************************
% subject optimal frequency selection based on AR FDA-F-Score
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
ORDER = 12; % order of AR model
NFFT = 1000;

muStep = 1; betaStep = 1;
muRange = [4,14]; betaRange = [16,40];
muSize = [4,5,6]; betaSize = [10,12];

n = length(data);

% caculate PSD ( power spectral density)
for i = 1:n
    [PSD{i}, f]=pburg(data{i},ORDER,NFFT,fs);
end
%%
% mu
nF = 0; % flag of frequency band
for kk = 1:length(muSize)
    for n_mu = 0 : floor( (muRange(2)-muRange(1)-muSize(kk)) / muStep )
        
        % PSD value selection
        mu0 = [muRange(1)+n_mu*muStep, muRange(1) + muSize(kk) + n_mu*muStep]
        rank = [find(f >= mu0(1), 1, 'first' ), find(f >= mu0(2), 1, 'first' )];
     
        n1 = 0;
        n2 = 0;
        % caculate feature value
        for i = 1:n
            muPSD = PSD{i}( rank(1) : rank(2) );
            if(labels(i) == 1)
                n1 = n1 + 1;
                BP1.mu(n1, :) = log10( var( muPSD ));
            else
                n2 = n2 + 1;
                BP2.mu(n2, :) = log10( var( muPSD ));
            end
        end
        % mean value
        BP1.mean = mean(BP1.mu);
        BP2.mean = mean(BP2.mu);
        % variance
        S1.mu = var( BP1.mu );
        S2.mu = var( BP2.mu);
        % caculate F-score
        nF = nF + 1;
        Res.FscoreMu(nF) = sum((BP1.mean - BP2.mean).^2) / sum(S1.mu + S2.mu);
        Res.mu(nF,:) = mu0;
    end
end
%%
% beta
nF = 0; % flag of frequency band
for kk = 1:length(betaSize)
    for n_beta = 0 : floor( (betaRange(2)-betaRange(1)-betaSize(kk)) / betaStep )
        
       % PSD value selection
        beta0 = [betaRange(1)+n_beta*betaStep, betaRange(1) + betaSize(kk) + n_beta*betaStep]
        rank = [find(f >= beta0(1), 1, 'first' ), find(f >= beta0(2), 1, 'first' )];
        
        n1 = 0;
        n2 = 0;
        % caculate feature value
        for i = 1:n
            betaPSD= PSD{i}( rank(1) : rank(2) );
            if(labels(i) == 1)
                n1 = n1 + 1;
                BP1.beta(n1, :) = log10( var( betaPSD ));
            else
                n2 = n2 + 1;
                BP2.beta(n2, :) = log10( var( betaPSD ) );
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