
function Rsquared = my_Rsquared_coeff(data,data_fit)
% R2 correlation coefficient computation

% The total sum of squares
sum_of_squares = sum((data-mean(data)).^2);

% The sum of squares of residuals, also called the residual sum of squares:
sum_of_squares_of_residuals = sum((data-data_fit).^2);

% definition of the coefficient of correlation is
Rsquared = 1 - sum_of_squares_of_residuals/sum_of_squares;

end
