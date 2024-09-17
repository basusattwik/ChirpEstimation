function obj = runLmcCore(obj)
%RUNLMCCORE Summary of this method goes here
%   Detailed explanation goes here

try 

    % Create a wait bar display
    waitbarIterCount = 1;
    wbar = waitbar(0, 'Sit tight!', 'Name','Running Langevin Monte Carlo', ...
                      'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    totalNumIter = obj.numIterNoise * obj.numIterLmc * obj.numParticles;
    
    % Initialize a live plotter
    if obj.bDisplayPlots
        ax = initLivePlots(obj);
    end

    totalGrads = zeros(obj.numParams, 1);
    for nind = 1:obj.numIterNoise

        % Save stepsize numerator
        stepSizeNum = obj.stepSizeInit .* (obj.noiseVar(1, nind) / obj.noiseVarFinal)^2;

        for tind = 1:obj.numIterLmc
            for pind = 1:obj.numParticles 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                                      %
                % --- Run MLE Model for Chirp Parameter Estimation --- %
                %                                                      %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Gradient smoothing
                for sind = 1:obj.numIterSmooth
                    extraNoise = obj.bEnableGaussSmooth * obj.noiseVar(1,nind) .* randn(obj.numParams,1);
                    paramNoisy = obj.param(:,pind) + extraNoise;
                    obj.cpe{1,pind} = obj.cpe{1,pind}.runCpeCore(paramNoisy); % This gives the gradients wrt all params for all particles
                    totalGrads = totalGrads + obj.cpe{1,pind}.dJ_phi.'; 
                end
                obj.grads(:,pind) = (totalGrads / obj.numIterSmooth);
                totalGrads(:) = 0;

                % Objective function for current param for current particle
                obj.cpe{1,pind}.J = obj.cpe{1,pind}.evalObjectiveFunc(obj.param(:,pind));  
                if obj.cpe{1,pind}.J <= obj.cpe{1,pind}.minObjTol
                    obj.cpe{1,pind}.bMinFound = true; % Check if J is small enough... that means we have converged
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                            %
                % --- Stepsize Annealing --- %
                %                            %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Compute gradient exp averaging         
                obj.avgGrads(:,pind) = obj.avgConst .* obj.grads(:,pind) + (1 - obj.avgConst) .* obj.avgGrads(:,pind);
                
                % Keep track of average gradient norm
                obj.gradNorm(1,pind)    = norm(obj.grads(:,pind))^2;
                obj.avgGradNorm(1,pind) = obj.avgConst .* abs(obj.gradNorm(1,pind)) + (1 - obj.avgConst) .* obj.avgGradNorm(1,pind);
                
                % Update stepsize
                obj.stepSize(:,pind) = stepSizeNum(:,pind) ./ (1 + obj.stepSizeConst .* obj.avgGrads(:,pind).^2);
                
                % Check if stepsize exceeds limits
                for npr = 1:obj.numParams
                    if obj.stepSize(npr,pind) > obj.stepSizeMax(npr,pind)
                        obj.stepSize(npr,pind) = obj.stepSizeMax(npr,pind);
                    end
                end    

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                             %
                % --- Simulated Tempering --- %
                %                             %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Adjust the temperature
                obj.temper = 1; %obj.tempConst * log10(1 + (tind));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                              %
                % --- Langevin Monte Carlo --- %
                %                              %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             

                % Do Langevin updates on all params to get a new proposed point
                paramProp = obj.param(:,pind) - obj.stepSize(:,pind) .* obj.grads(:,pind) + ...
                                                obj.bEnableLangevin * sqrt(2 * obj.stepSize(:,pind) ./ obj.temper) .* randn(obj.numParams,1);

                obj.saveProposedParams(:,pind, tind, nind) = paramProp;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                   %
                % --- Run Metropolis Correction --- %
                %                                   %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % obj.param(:,pind) = paramProp;
                alpha = min(1, exp(-obj.temper * (obj.cpe{1,pind}.evalObjectiveFunc(paramProp) - obj.cpe{1,pind}.J)));
                if rand(1,1) <= alpha
                    obj.param(:,pind) = paramProp;
                    bAccept = true;
                else
                    bAccept = false;
                end

                % Add proposal density


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                %
                % --- Save Data for Analysis --- %
                %                                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
                % Fill analysis arrays
                obj.saveParams(:,pind, tind, nind)    = obj.param(:,pind);
                obj.saveGrads(:, pind, tind, nind)    = obj.grads(:,pind);
                obj.saveStepSize(:,pind,tind,nind)    = obj.stepSize(:,pind); 
                obj.saveGradNorm(pind, tind, nind)    = obj.gradNorm(1,pind);
                obj.saveAvgGradNorm(pind, tind, nind) = obj.avgGradNorm(1,pind);
                obj.saveObjFunc(pind, tind, nind)     = obj.cpe{1,pind}.J;
                obj.savebAccept(pind,tind,nind)       = bAccept;

                %%%%%%%%%%%%%%%%%%%
                %                 %
                % --- Waitbar --- %
                %                 %
                %%%%%%%%%%%%%%%%%%%

                % Increment iteration count for waitbar
                waitbarIterCount = waitbarIterCount + 1;

                % Update waitbar and message
                if getappdata(wbar, 'canceling')
                    disp('Simulation cancelled!')
                    obj.bStopSim = true;
                    delete(wbar);
                    return
                end    
                if mod(waitbarIterCount, 10) == 0
                    waitbar(waitbarIterCount / totalNumIter, wbar, ['Noise Iteration: ', num2str(nind), '/', num2str(obj.numIterNoise)]);
                end

            end % end numParticles
          
            % Save temperature
            obj.saveTemper(tind,nind) = obj.temper;

            % Stopping criterion
            for pind = 1:obj.numParticles
                if obj.cpe{1,pind}.bMinFound
                    obj.bStopLmc = true;
                    break; 
                end
            end
        end % end numIterLmc

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                     %
        % --- Starting point for the next noise iteration --- %
        %                                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        % Get the time index for which the param led to the lowest objective func value
        for pind = 1:obj.numParticles
            [~, minObjLmcInd] = min(obj.saveObjFunc(pind, 1:tind, nind));
            obj.param(:,pind) = obj.saveParams(:,pind, minObjLmcInd, nind); % This is the starting point for the next round of LMC updates
        end

        % Update Live Plot
        if obj.bDisplayPlots
            updateLivePlots(obj, tind, nind, ax);
        end

        % Stop iterations since min has been found
        if obj.bStopLmc
            disp(['Minimum has been found already. LMC stopped after ', num2str(nind), ' iterations.']);
            break;
        end
        
    end % end numIterNoise

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                           %
    % --- Find the optimum parameter vector --- %
    %                                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    objFuncVals = zeros(obj.numParticles, 1);
    for pind = 1:obj.numParticles
        objFuncVals(pind, 1) = obj.cpe{1,pind}.evalObjectiveFunc(obj.param(:,pind));
    end

    [~, obj.bestParticleInd] = min(objFuncVals);
    obj.optParam = obj.param(:, obj.bestParticleInd);

    % Finally, compute the scalar gains for CPE
    obj.cpe{1,obj.bestParticleInd} = obj.cpe{1,obj.bestParticleInd}.compScalarGains(obj.optParam);   

    % Reconstruct the chirp from the estimated parameters
    obj.cpe{1,obj.bestParticleInd} = obj.cpe{1,obj.bestParticleInd}.reconChirpSignals();

    % Evaluate errors in estimation
    obj.cpe{1,obj.bestParticleInd} = obj.cpe{1,obj.bestParticleInd}.evalParamErrors();

    % Housekeeping
    disp('Simulation complete!')
    delete(wbar);

    fprintf('\n');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['The estimated phase parameters are     =  ', num2str(obj.optParam.')]);
    disp(['The estimated amplitude parameters are =  ', num2str(obj.cpe{1,obj.bestParticleInd}.rhoEst')]);
    fprintf('\n');
    disp(['Found by particle ', num2str(obj.bestParticleInd)]);
    fprintf('\n');
    disp(['Log Error for phase params     = ', num2str(log10(obj.cpe{1,obj.bestParticleInd}.sqrPhiError.'))]);
    disp(['Log Error for amplitude params = ', num2str(log10(obj.cpe{1,obj.bestParticleInd}.sqrRhoError.'))]);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf('\n');

    % Housekeeping
    disp('Simulation complete!')
    delete(wbar);

catch me
    delete(wbar); % Close wait bar when simulation errors out
    throw(me);
end

end