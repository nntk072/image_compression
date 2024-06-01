% Reconstruction function
function Brec = reconstruct_image(E)
    Brec = zeros(size(E));
    for i = 1:size(E, 1)
        for j = 1:size(E, 2)
            predictors = [];
            if i > 1
                predictors = [predictors, Brec(i-1, j)];
            end
            if j > 1
                predictors = [predictors, Brec(i, j-1)];
            end
            if i > 1 && j > 1
                predictors = [predictors, Brec(i-1, j-1)];
            end
            if isempty(predictors)
                y_MAP = 128;
            else
                y_MAP = median(predictors);
            end
            Brec(i, j) = E(i, j) + round(y_MAP);
        end
    end
end