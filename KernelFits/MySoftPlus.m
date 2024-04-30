function [DataOut] = MySoftPlus(DataIn, whichDim)

if nargin<2
    whichDim = 1;
end

DataOut = log(1 + exp(DataIn));

end