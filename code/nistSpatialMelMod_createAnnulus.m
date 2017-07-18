function [ weight ] = nistSpatialMelMod_createAnnulus( eccentricityDeg, radiusInnerEdgeAnnulusDeg, radiusOuterEdgeAnnulusDeg, widthHalfCosineSmoothDeg )
% function [ weight ] = nistSpatialMelMod_createAnnulus( eccentricityDeg, radiusInnerEdgeAnnulusDeg, radiusOuterEdgeAnnulusDeg, widthHalfCosineSmoothDeg )
%
% Implements an annulus with half-cosine edges, under the control of the
% passed parameters. The output is the weight (from 0 to 1) to be applied
% to stimulus contrast at a point on the display corresponding to the
% passed eccentricityDeg.
%

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