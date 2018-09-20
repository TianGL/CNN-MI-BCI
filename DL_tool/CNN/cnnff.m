function net = cnnff(net, x)
    n = numel(net.layers);
    net.layers{1}.a{1} = x;
    inputmaps = 1;

    for r = 2 : n   %  for each layer
        if strcmp(net.layers{r}.type, 'c')
            %  !!below can probably be handled by insane matrix operations
            for j = 1 : net.layers{r}.outputmaps   %  for each output map
                %  create temp output map
                z = zeros(size(net.layers{r - 1}.a{1}) - [net.layers{r}.kernelsize(1) - 1 net.layers{r}.kernelsize(2) - 1 0]);
                for i = 1 : inputmaps    %  for each input map
                    %  convolve with corresponding kernel and add to temp output map
                    z = z + convn(net.layers{r - 1}.a{i}, net.layers{r}.k{i}{j}, 'valid');
                end
                %  add bias, pass through nonlinearity 
                net.layers{r}.a{j} = ReLU(z + net.layers{r}.b{j});% ReLU function
%                 net.layers{r}.a{j} = sigm(z + net.layers{r}.b{j}); % sigmod function
            end
            %  set number of input maps to this layers number of outputmaps
            inputmaps = net.layers{r}.outputmaps;
        elseif strcmp(net.layers{r}.type, 's')
            %  downsample
            for j = 1 : inputmaps
%                 z = convn(net.layers{r - 1}.a{j}, ones(1, net.layers{r}.scale,1) / (net.layers{r}.scale ^ 2), 'valid');   %  !! replace with variable
%                 net.layers{r}.a{j} = z(1 : net.layers{r}.scale : end, 1 : net.layers{r}.scale : end, :);
                 [z , pos] = maxpooling(net.layers{r-1}.a{j}, net.layers{r}.scale);
                 net.layers{r}.a{j} = z;
                 net.layers{r}.maxpos{j} = pos;
            end
        end
    end

    %  concatenate all end layer feature maps into vector
    net.fv = [];
    for j = 1 : numel(net.layers{n}.a)
        sa = size(net.layers{n}.a{j});
        net.fv = [net.fv; reshape(net.layers{n}.a{j}, sa(1) * sa(2), sa(3))];
    end
    %  feedforward into output perceptrons
    net.o = ReLU(net.ffW * net.fv + repmat(net.ffb, 1, size(net.fv, 2))); % ReLU
%     net.o = sigm(net.ffW * net.fv + repmat(net.ffb, 1, size(net.fv, 2))); % sigmod
end
