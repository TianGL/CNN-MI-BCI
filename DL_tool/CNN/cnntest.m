function [er, bad, predict, origin] = cnntest(net, x, y)
    %  feedforward
    net = cnnff(net, x);
    [~, predict] = max(net.o);
    [~, origin] = max(y);
    bad = find(predict ~= origin);

    er = numel(bad) / size(y, 2);
end
