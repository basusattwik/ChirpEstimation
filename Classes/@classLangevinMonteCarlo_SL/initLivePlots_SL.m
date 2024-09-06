function ax = initLivePlots_SL(obj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% params = [];
% for c = 1:obj.cpe{1,1}.Nc
%     params = [params; obj.cpe{1,1}.phi{1,c}(2:end,1)]; %#ok<AGROW>
% end
% minObjFuncVal = obj.cpe{1,1}.evalObjectiveFunc(params);

figure('windowstyle','docked');
    ax = plot(1:obj.numIterLmc, -ones(obj.numIterLmc,obj.numParticles),  'LineWidth', 1.1); hold on;
    % plot(1:obj.numIterLmc, minObjFuncVal * ones(obj.numIterLmc,1), '-.', 'LineWidth', 1.1, 'DisplayName', 'Min Objective Func.');
    grid on; grid minor; 
    xlabel('Iterations', 'FontSize', 12); 
    ylabel('Objective Func', 'FontSize', 12); 
    xlim([0, obj.numIterLmc]);
    ylim([0, inf]);
    title('Objective Function vs Iterations', 'FontSize', 14);
    
end