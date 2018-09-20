function [er,predict, origin] = trailtest(net, data, testflag)

    % get data for each trial
    n = 1;
    trailS = length(data.tra{1}.image);
    for i = 1:length(data.tra)
        if(testflag(i))
            for j = 1:trailS
                testData{n}(:,:,j) = data.tra{i}.image{j};
                testLabels{n}(:,j) =  data.tra{i}.labels{j};
            end
            n = n+1;
        end
    end
    
    right = 0;
    for i = 1:n-1
        x = testData{i};
        y = testLabels{i};
        %  feedforward
        net = cnnff(net, x);
        [~, pre] = max(net.o);
        [~, ori] = max(y);
        origin(i) = ori(1,1,1);
        wrong = find(pre ~= ori);
        traER = numel(wrong) / size(y, 2);
        if(traER<0.5)
            right = right+1;
            predict(i) = origin(i);
        else
            if(origin(i) == 1)
                predict(i) = 2;
            else
                predict(i) = 1;
            end
        end
    end
    er = 1 - right/length(testData);
end
