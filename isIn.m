function [ hit ] = isIn( INDEXG,T )
%ISIN Summary of this function goes here
%   Detailed explanation goes here
hit=false;
for i=1:size(INDEXG,2)
    if T==INDEXG(i)
        hit=true;
    end
end
end

