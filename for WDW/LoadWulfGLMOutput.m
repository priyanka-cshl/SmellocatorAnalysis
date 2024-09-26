function [PSTHOut] = LoadWulfGLMOutput(myPath,nUnits)
load(myPath);

for i = 1:nUnits
    PSTHOut(i,:) = eval(['neuron_',num2str(i-1),'.full_test_y_pred']);
end

end