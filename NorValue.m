function res = NorValue(array, ptg)
% normalization input data
% array(): original data
% ptg: percentage of max value of normalization data in orginal data

    temp = array(:);
    temp = sort(temp);
    maxValue = temp(floor(length(temp) * ptg));
    deNumber = find(array>maxValue);
    array(deNumber) = maxValue;
    res = array ./ maxValue;
end