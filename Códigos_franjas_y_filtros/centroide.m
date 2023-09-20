
function [R, C] = centroide(I)
C=round((1:size(I,2))*sum(I,1)'/sum(I(:)));
R=round((1:size(I,1))*sum(I,2)/sum(I(:)));

