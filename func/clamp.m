function [in] = clamp(in,low,high)
%clamp limit a certain input to given bounds
%   Detailed explanation goes here
if in > high
        in = high;
    elseif in < low
        in = low;
end
end