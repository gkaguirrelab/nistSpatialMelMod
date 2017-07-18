function [ weight ] = nistSpatialMelMod_createAnnulus( eccentricityDeg, radiusInnerEdgeAnnulusDeg, radiusOuterEdgeAnnulusDeg, widthHalfCosineSmoothDeg )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if eccentricityDeg < radiusInnerEdgeAnnulusDeg
    weight = 0;
    return
end

if eccentricityDeg > radiusOuterEdgeAnnulusDeg
    weight = 0;
    return
end

if (eccentricityDeg > radiusInnerEdgeAnnulusDeg+widthHalfCosineSmoothDeg) && ...
        (eccentricityDeg < radiusOuterEdgeAnnulusDeg-widthHalfCosineSmoothDeg)
    weight = 1;
    return
end

if (eccentricityDeg >= radiusInnerEdgeAnnulusDeg) && ...
        (eccentricityDeg <= radiusInnerEdgeAnnulusDeg+widthHalfCosineSmoothDeg)
    weight = cos( ((eccentricityDeg-radiusInnerEdgeAnnulusDeg)/widthHalfCosineSmoothDeg) * pi + pi);
    weight = (weight + 1)/2;
    return
end

if (eccentricityDeg >= radiusOuterEdgeAnnulusDeg-widthHalfCosineSmoothDeg) && ...
        (eccentricityDeg <= radiusOuterEdgeAnnulusDeg)
    weight = cos( ((eccentricityDeg-radiusOuterEdgeAnnulusDeg-widthHalfCosineSmoothDeg)/widthHalfCosineSmoothDeg) * pi );
    weight = (weight + 1)/2;
    return
end

error('The eccentricity did not meet any conditions for the mask');

end % function