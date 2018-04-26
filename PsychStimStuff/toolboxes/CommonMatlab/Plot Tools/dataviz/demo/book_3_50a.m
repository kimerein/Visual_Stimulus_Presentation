%  book_3_50.m
%  calls normalqqplot

load fly

uniqueTemp = unique(Temperature);
nT = length(uniqueTemp);

for ii = 1:nT
   index = Temperature==uniqueTemp(ii);
   meanFacet(ii) = mean(FacetNumber(index));
end

p = polyfit(uniqueTemp(:),meanFacet(:),1);
fitFacet = polyval(p,uniqueTemp);
residual = meanFacet(:)-fitFacet(:);

normalqqplot(residual)

return
xlabel('Temperature (deg C)')
ylabel('Residual Facet Number')
title('Fly')
