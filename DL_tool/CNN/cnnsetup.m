function net = cnnsetup(net, x, y)
%     assert(~isOctave() || compare_versions(OCTAVE_VERSION, '3.8.0', '>='), ['Octave 3.8.0 or greater is required for CNNs as there is a bug in convolution in previous versions. See http://savannah.gnu.org/bugs/?39314. Your version is ' myOctaveVersion]);
    inputmaps = 1;
    mapsize = size(squeeze(x(:, :, 1))); % get size of input image

    for r = 1 : numel(net.layers)   %  layer
        if strcmp(net.layers{r}.type, 's')
            mapsize = mapsize ./ [1,net.layers{r}.scale];
            assert(all(floor(mapsize)==mapsize), ['Layer ' num2str(r) ' size must be integer. Actual: ' num2str(mapsize)]);
            for j = 1 : inputmaps
                net.layers{r}.b{j} = 0;
            end
        end
        if strcmp(net.layers{r}.type, 'c')
            mapsize = mapsize - net.layers{r}.kernelsize + 1;
            fan_out = net.layers{r}.outputmaps * net.layers{r}.kernelsize(1) * net.layers{r}.kernelsize(2);
            for j = 1 : net.layers{r}.outputmaps  %  output map
                fan_in = inputmaps * net.layers{r}.kernelsize(1) * net.layers{r}.kernelsize(2);
                for i = 1 : inputmaps  %  input map
                    net.layers{r}.k{i}{j} = (rand(net.layers{r}.kernelsize) - 0.5) * 2 * sqrt(6 / (fan_in + fan_out));
                end
                net.layers{r}.b{j} = 0;
            end
            inputmaps = net.layers{r}.outputmaps;
        end
    end
    % 'onum' is the number of labels, that's why it is calculated using size(y, 1). If you have 20 labels so the output of the network will be 20 neurons.
    % 'fvnum' is the number of output neurons at the last layer, the layer just before the output layer.
    % 'ffb' is the biases of the output neurons.
    % 'ffW' is the weights between the last layer and the output neurons. Note that the last layer is fully connected to the output layer, that's why the size of the weights is (onum * fvnum)
    fvnum = prod(mapsize) * inputmaps; %number of output neuron
    onum = size(y, 1); % output

    net.ffb = zeros(onum, 1);
    net.ffW = (rand(onum, fvnum) - 0.5) * 2 * sqrt(6 / (onum + fvnum));
end
