function bounds = conf(data,low,high)
%low and high should be given in decimal values not percents

sd = sort(data);
len = length(sd);

low_idx = floor(len*low);
high_idx = ceil(len*high);

bounds = [sd(low_idx),sd(high_idx)];
end