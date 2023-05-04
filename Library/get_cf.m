function [CF_List] = get_cf()
%get_cf() prompts the user to enter CF factors
%   CF factors should be normalized and CF1 should always be equal to one. 
fprintf('\nEnter CF factors\n')
fprintf('CF1==> 1 ');
CF1 = 1;
CF2 = input('\nCF2==>  ');
CF3 = input('\nCF3==>  ');
CF4 = input('\nCF4==>  ');
CF_List = [CF1 CF2 CF3 CF4];
end

