function [CantMeas,Cant2002]=SimpleCantEstimateLR(Longitude,Latitude,Depth,Dates);
    % Uses output from Siv Lauvset et al's gridded 2002 TTD-based Cant
    % estimates for 2002 to grid to desired locations;
    % DO NOT USE THIS TO ESTIMATE CANTH FOR ANY OTHER APPLICATIONS.
    % There is undoubtedly a better way.
    tic
    load Cant2002LEAInterpolantLR.mat  
    % We want latitude distances to count for more than longitude distances
    % and we need to account for differing vertical units (meters) and
    % horizontal units (degrees).  The ScalingConstants help with that.
    InterpolantLR.Points(:,1)=InterpolantLR.Points(:,1)*InterpScalingConsts.Longitude;
    InterpolantLR.Points(:,3)=InterpolantLR.Points(:,3)*InterpScalingConsts.Depth;
    Cant2002=InterpolantLR(Longitude*InterpScalingConsts.Longitude,Latitude,Depth*InterpScalingConsts.Depth);
    % Here we borrow from Gruber et al. 2019's analysis to simplistically
    % scale Cant2002 to an arbitrary year by assuming transient steady
    % state and an exponential increase in pCO2 based on the 1994 to 2007
    % increase.  There would be a better way of doing this calculation if
    % we knew what the future pCO2 would be... but if we knew that we'd
    % have other papers to write.  This could be improved using a full pCO2
    % historical record, but this simplistic assumption is already a
    % meaningful improvement over alternative simplistic approaches for
    % dealing with OA (e.g. LIPHR's).
    CantMeas=Cant2002.*exp(0.018989*(Dates-2002));
end