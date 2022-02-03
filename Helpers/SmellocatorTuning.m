function [myCurve,XBins] = SmellocatorTuning(VarType, VarX, VarY)

switch VarType
    case 'Lever'
        Binsize = 0.25;
        XBins(:,1) = 0:Binsize:(5-Binsize);
        XBins(:,2) = XBins(:,1) + Binsize;
    case 'Odor'
        Binsize = 10;
        XBins(:,1) = 20:Binsize:(230-Binsize);
        XBins(:,2) = XBins(:,1) + Binsize;
end

% VarX can be a column vector, or a matrix
% VarY can also be a column vector, or a matrix
myCurve = NaN*zeros(size(XBins,1),4,size(VarX,2),size(VarY,2));
for whichUnit = 1:size(VarY,2)
    YVar = VarY(:,whichUnit);
    for whichVar = 1:size(VarX,2)
        XVar = VarX(:,whichVar);
        for myBin = 1:size(XBins,1)
            idx = find((XVar>=XBins(myBin,1))&(XVar<XBins(myBin,2)));
            if ~isempty(idx)
                myCurve(myBin,:,whichVar,whichUnit) = [mean(YVar(idx)) median(YVar(idx)) std(YVar(idx))/sqrt(numel(idx)) numel(idx)];
            end
        end
    end
end
end

