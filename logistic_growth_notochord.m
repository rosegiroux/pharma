% Rose Long
% Feb 2 2024
% use logistic model and measured mitotic index
%  to predict polulation at later timepoints
% (dy/dx) = r * y * (1 - (y/K))
% 
% where r is the growth rate and K is the carrying capacity.


%Define the function for the logistic growth model:
function dy = logistic_growth_notochord(t, y, r, K)
    dy = r * y * (1 - y/K);
end

