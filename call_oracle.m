function [nextseed, y, feasibilityFLAG]= call_oracle(nextseed, x, feascheckonlyFLAG,OracleName)
    if strcmp(OracleName,'BusScheduling') == 1
        [nextseed,y,feasibilityFLAG] = OracleBusWaitingTime(nextseed,x,feascheckonlyFLAG);
    elseif strcmp(OracleName,'DiscreteQuadratic') == 1
        [nextseed,y,feasibilityFLAG] = OracleDiscreteQuadratic(nextseed,x,feascheckonlyFLAG);
    elseif strcmp(OracleName,'DynamicNews') == 1
        [nextseed,y,feasibilityFLAG] = OracleDynamicNews(nextseed,x,feascheckonlyFLAG);
    end
end