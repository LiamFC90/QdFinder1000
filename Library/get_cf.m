function [CF_List] = get_cf()
%get_cf() prompts the user to enter CF factors
%   CF factors should be normalized and CF1 should always be equal to one. 
fprintf('\nEnter CF factors\n')
fprintf('CF1==> 1.00\n');
CF1 = 1;
CF2 = input('CF2==>  .')/100;
CF3 = input('CF3==>  .')/100;
CF4 = input('CF4==>  .')/100;
CF_List = [CF1 CF2 CF3 CF4];
end

