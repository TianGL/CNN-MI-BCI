% prerequist running ~/biosig4octmat-3.2.0/biosig_installer.m

%{
 @ name: Classificaiton of left-right hands motor imagery eeg signals
 @ author: Geliang Tian
 @ email: tglasd@163.com
 @ year: 2018
 @ version: 1.0
%}

%********************************************************************************
% 1) Read the data from dataset and get the target segments data
clear all
clc

% dictionary of dataset and label
addpath('D:\BCI\EEG dataset\bci competition IV 2b\dataset\');
addpath('D:\BCI\EEG dataset\bci competition IV 2b\dataset\true_labels');

% get data of each subject
for kn=1:9
    p(kn) = 0; % get the number of none value
    sampleRate = 250; % sample frequency
    
    % to enrich data, set window size to 2 s with an overlap of 90%
    delayTime = 50;
    timeScale = 2;
    
    total_data_idx = 1; % combinitation of data from all session; for optimal frequency band selection
    for i = 1:3 % set 5 to also get evalution dataset
        if(i<4)
            dataName = ['B0',num2str(kn),'0',num2str(i),'T.gdf'];
            labelName = ['B0',num2str(kn),'0',num2str(i),'T.mat'];
        else
            dataName = ['B0',num2str(kn),'0',num2str(i),'E.gdf'];
            labelName = ['B0',num2str(kn),'0',num2str(i),'E.mat'];
        end
        [signal{i},H{i}] = sload(dataName);
        load(labelName);
        tureLabels{i}= classlabel;
        
        % get the label of each segments
        CIV2b_S{kn}.D{i}.raw = signal{i}(:,1:3); % original signals from c3,c4 and cz
        trail = length(H{i}.TRIG);
        
        CIV2b_S{kn}.D{i}.labels = tureLabels{i}; % get label
        
        n = 1;
        for j = 1:trail % extract MI signals from 3-7s
            meanValue = mean( CIV2b_S{kn}.D{i}.raw(H{i}.TRIG(j) : H{i}.TRIG(j) + sampleRate*2 ,:));
            if( any( any( isnan( CIV2b_S{kn}.D{i}.raw( H{i}.TRIG(j) : H{i}.TRIG(j) + sampleRate*7 ,: ) ) ) ) )
                p(kn) = p(kn) + 1;
                continue;
            end
            temp = CIV2b_S{kn}.D{i}.raw( H{i}.TRIG(j) + sampleRate*3 : H{i}.TRIG(j) + sampleRate*7 ,: );
            CIV2b_S{kn}.D{i}.MI{j}= [temp(:,1) - meanValue(1),temp(:,2) - meanValue(2),temp(:,3) - meanValue(3)];
            m = 0;
            while ( m*delayTime <= 500) % get the 2s length segment
                CIV2b_S{kn}.D{i}.tra{n}.data{m+1} =  CIV2b_S{kn}.D{i}.MI{j}(m * delayTime + 1 : m * delayTime + timeScale * sampleRate,:);
                CIV2b_S{kn}.D{i}.tra{n}.labels{m+1} =  CIV2b_S{kn}.D{i}.labels(j);
                
                % data combination
                CIV2b_S{kn}.all_data{total_data_idx} = CIV2b_S{kn}.D{i}.tra{n}.data{m+1};
                CIV2b_S{kn}.all_data_label(total_data_idx) = CIV2b_S{kn}.D{i}.tra{n}.labels{m+1};
                total_data_idx = total_data_idx + 1;
                
                m = m+1;
            end
            CIV2b_S{kn}.D{i}.Labels(n) = CIV2b_S{kn}.D{i}.labels(j);
            n = n+1;
        end
    end
end

clearvars -EXCEPT CIV2b_S p
%%
%**********************************************************
% 2) preprocessing and form of input image of CNN
for kn = 1:9
    
    %subject optimal frequency selection BP FDA-F-Score
%     [muBand, betaBand] = BPFeatureBandSelection(CIV2b_S{kn}.all_data, CIV2b_S{kn}.all_data_label, 250);
    
    % subject optimal frequency selection AR PSD FDA-F-Score
%     [muBand, betaBand] = ARFeatureBandSelection(CIV2b_S{kn}.all_data, CIV2b_S{kn}.all_data_label, 250);
    
    for sess = 1:3 % 5 for contains evaluation session
        % get training data and test data
        CIV2b_Data_S{kn}.se{sess}.Labels = CIV2b_S{kn}.D{sess}.Labels;
        
        % extend frequency bands (better performance)
        muBand = [4,14];
        betaBand = [16,32];
        
        % get the all input images and labelsof CNN
        for i = 1 : length(CIV2b_S{kn}.D{sess}.tra)
            for kk = 1:length(CIV2b_S{kn}.D{sess}.tra{i}.data)
                for j = 1:3
                    Cx{j} = CIV2b_S{kn}.D{sess}.tra{i}.data{kk}(:,j);
                    
                    % short time Fourier transform
                    [Fstft, f, t] = stft(Cx{j}, 64, 14, 512, 250);
                    Mu{j} = abs( Fstft( (find(f<muBand(1),1,'last') ) : (find(f<muBand(2),1,'last')) +1 ,:) );
                    Beta = abs( Fstft( (find(f<betaBand(1),1,'last') ) : (find(f<betaBand(2),1,'last')) +1 ) );%%和paper不一样
                    
                    % beta band cubic interpolation
                    interNum = size(Mu{j},1);
                    fBeta = betaBand(1) : (betaBand(2)-betaBand(1))/(interNum-1) : betaBand(2);
                    [X,Y] = meshgrid(t,f);
                    [X1,Y1] = meshgrid(t,fBeta);
                    Beta_intrp{j} = interp2( X,Y,abs( Fstft ),X1,Y1,'cubic');
                    
                    % normalization
                    Mu{j} = NorValue(Mu{j},1);
                    Beta_intrp{j} = NorValue(Beta_intrp{j}, 1);
                end
                
                CIV2b_Data_S{kn}.se{sess}.tra{i}.C3{kk} = [Beta_intrp{1}; Mu{1}];
                CIV2b_Data_S{kn}.se{sess}.tra{i}.Cz{kk} = [Beta_intrp{2}; Mu{2}];
                CIV2b_Data_S{kn}.se{sess}.tra{i}.C4{kk} = [Beta_intrp{3}; Mu{3}];
                
                CIV2b_Data_S{kn}.se{sess}.tra{i}.image{kk} =  [CIV2b_Data_S{kn}.se{sess}.tra{i}.C4{kk};CIV2b_Data_S{kn}.se{sess}.tra{i}.Cz{kk};CIV2b_Data_S{kn}.se{sess}.tra{i}.C3{kk}];
                switch CIV2b_S{kn}.D{sess}.tra{i}.labels{kk} % for each label
                    case 1
                        CIV2b_Data_S{kn}.se{sess}.tra{i}.labels{kk} = [1; 0];
                    case 2
                        CIV2b_Data_S{kn}.se{sess}.tra{i}.labels{kk} = [0; 1];
                    otherwise
                        CIV2b_Data_S{kn}.se{sess}.tra{i}.labels{kk} = [0; 0];
                end
            end
        end
        
        CIV2b_Data_S{kn}.band = [muBand,betaBand];
    end
end

clearvars -EXCEPT CIV2b_Data_S

%%
%************************************************
% 3) Convolution Neuro Network training and evaluation
%
% Deep learning toolbox contributor information
% http://www2.imm.dtu.dk/pubdb/views/publication_details.php?id=6284
% @MASTERSTHESIS\{IMM2012-06284,
%     author       = "R. B. Palm",
%     title        = "Prediction as a candidate for learning deep hierarchical models of data",
%     year         = "2012",
% }
% Contact: rasmusbergpalm at gmail dot com

% add path of DL code
addpath(genpath('D:\BCI\BCI code\program1\matlab\DL_tool'))

for kn = 1:9
    disp(['CNN_kn=' num2str(kn)]);
    % network structure
    kernlSize = [ size(CIV2b_Data_S{kn}.se{1}.tra{1}.image{1}, 1), 3];
    cnn.layers = {
        struct('type', 'i') %input layer
        struct('type', 'c', 'outputmaps', 30, 'kernelsize', kernlSize) %convolution layer
        struct('type', 's', 'scale', 10) %sub sampling layer
        };
    
    % initialization parameters of CNN: 0.02 40 1500
    opts.alpha = 0.05;
    opts.batchsize = 40; %40,50
    opts.numepochs = 300;
    
    indices = [];
    dataLabels = [];
    data = [];
    for sess = 1:3 % 5
        % get training dataset and test dataset
        sess_datalen = length(CIV2b_Data_S{kn}.se{sess}.Labels);
        sess_dataLabels = zeros(sess_datalen, 2);
        for i = 1:sess_datalen
            switch CIV2b_Data_S{kn}.se{sess}.Labels(i)
                case 1
                    sess_dataLabels(i,:) = [1,0];
                case 2
                    sess_dataLabels(i,:) = [0,1];
                otherwise
                    sess_dataLabels(i,:) = [0,0];
            end
        end
        
        %10-fold cross-validation
        indices_s= crossvalind('Kfold',sess_dataLabels(:,1), 10);
        
        lenCur = length(indices);
        indices( lenCur+1 : lenCur+sess_datalen, 1 ) = indices_s;
        dataLabels( lenCur+1 : lenCur+sess_datalen, : ) = sess_dataLabels;
        data.tra( lenCur+1 : lenCur+sess_datalen) = CIV2b_Data_S{kn}.se{sess}.tra;
        data.Labels( lenCur+1 : lenCur+sess_datalen) = CIV2b_Data_S{kn}.se{sess}.Labels;
    end
    
    % training and test each fold
    for fold = 1:10
        test = (indices == fold); train = ~test;
        lenTrain = length( find(train == 1)) * length(data.tra{1}.image);
        lenTest = length( find(test == 1)) * length(data.tra{1}.image);
        
        train_data = zeros(kernlSize(1),32,lenTrain);
        train_labels = zeros(2,lenTrain);
        test_data = zeros(kernlSize(1),32,lenTest);
        test_labels = zeros(2,lenTest);
        
        trailS = length(data.tra{1}.image); % size of sample of each trial
        
        % combine the training data from each session
        n = 1;
        m = 1;
        for i = 1:sess_datalen
            if(train(i) == 1)
                for j = 1:trailS
                    train_data(:,:,n) = data.tra{i}.image{j};
                    train_labels(:,n) = data.tra{i}.labels{j};
                    n = n+1;
                end
            else
                for j = 1:trailS
                    test_data(:,:,m) = data.tra{i}.image{j};
                    test_labels(:,m) = data.tra{i}.labels{j};
                    m = m+1;
                end
            end
        end
        
        % regard for the batch size, remove little data from training set
        len = floor(size(train_data,3)/opts.batchsize);
        train_data = train_data(:,:,1:len*opts.batchsize);
        train_labels = train_labels(:,1:len*opts.batchsize);
        
        %CNN setup and training
        cnn = cnnsetup(cnn, train_data, train_labels);
        cnn = cnntrain(cnn, train_data, train_labels, opts);
        
        [er, predict, origin] = trailtest(cnn, data, test); % eavuation
        
        % evalutaion result
        resfinal_s13{kn}.F{fold}.er = er;
        resfinal_s13{kn}.S{fold}.predict = predict;
        resfinal_s13{kn}.S{fold}.origin = origin;
        ERR(kn,fold) = er;
    end
end