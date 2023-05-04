function [] = posstats(posdata)

%This function finds and prints the spacing and associated error. 
%spacing then error. Spacing is found using sqrt(x^2 + y^2) and error is
%likewise found.
%

pos = sqrt(posdata(4,2)^2+posdata(4,5)^2);
error = sqrt(posdata(4,4)^2+posdata(4,7)^2);
fprintf('\nSpacing is %f\n',pos)
fprintf('error is %f\n',error)

end