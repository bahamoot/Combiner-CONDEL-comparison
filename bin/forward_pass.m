function [out1, out2, out3] = forward_pass(weight_level1, weight_level2, weight_level3, data)

datasize = size(data, 2);

in1  = weight_level1 * [data ; ones(1, datasize)];
out1 = [2 ./ (1+exp(-in1)) - 1 ; ones(1, datasize)];
in2  = weight_level2 * out1;
out2 = [2 ./ (1+exp(-in2)) - 1 ; ones(1, datasize)];
in3  = weight_level3 * out2;
out3 = 1 ./ (1+exp(-in3));
