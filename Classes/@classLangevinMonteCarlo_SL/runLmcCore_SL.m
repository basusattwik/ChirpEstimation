function obj = runLmcCore_SL(obj)
    %RUNLMCCORE Summary of this method goes here
    %   Detailed explanation goes here

    try 

        % Create a wait bar display
        wbar = waitbar(0, 'Sit tight!', 'Name','Running Langevin Monte Carlo', ...
                          'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
        
        % Initialize a live plotter
        if obj.bDisplayPlots
            ax = initLivePlots_SL(obj);
        end

        totalGrads = zeros(obj.numParams, 1);
        totalHess  = zeros(obj.numParams, obj.numParams);
        for tind = 1:obj.numIterLmc

            % Save temperature
            obj.saveTemper(tind,1) = obj.temper;

            for pind = 1:obj.numParticles 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                                      %
                % --- Run MLE Model for Chirp Parameter Estimation --- %
                %                                                      %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                obj.cpe{1,pind} = obj.cpe{1,pind}.runCpeCore(obj.param(:,pind));
                currentFuncVal  = obj.cpe{1,pind}.J; % Function value at the unperturbed parameter

                % Gradient smoothing and Hessian approximation
                for sind = 1:obj.numIterSmooth

                    normRand   = randn(obj.numParams,1);
                    extraNoise = obj.bEnableGaussSmooth * obj.noiseVar(1,pind) .* normRand;
                    paramNoisy = obj.param(:,pind) + extraNoise;

                    obj.cpe{1,pind} = obj.cpe{1,pind}.runCpeCore(paramNoisy); % This gives the gradients wrt all params for all particles
                    
                    % Keep adding gradients
                    totalGrads = totalGrads + obj.cpe{1,pind}.dJ_phi.'; 

                    % Approximate Hessian and keep adding. This uses
                    % Stein's lemma
                    totalHess  = totalHess  + (normRand * normRand' -  obj.Id) * (obj.cpe{1,pind}.J - currentFuncVal) * (1/obj.noiseVar(1,pind)^2);
                end
                obj.grads(:,pind)  = (totalGrads / obj.numIterSmooth);
                obj.hess(:,:,pind) = (totalHess  / obj.numIterSmooth);
                totalGrads(:) = 0;
                totalHess(:)  = 0;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                          %
                % --- Update Noise Var --- %
                %                          %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Fill in Hessian approx
                obj.absHessianTrace(1,pind) = obj.avgConst * abs(trace(obj.hess(:,:,pind))) + (1 - obj.avgConst) * obj.absHessianTrace(1,pind);
                obj.noiseVar(1,pind)        = min(obj.noiseVarInit, max(obj.noiseVarMin, obj.noiseVar(1,pind) - obj.stepNoiseVar(1,pind) * obj.absHessianTrace(1,pind)));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                            %
                % --- Stepsize Annealing --- %
                %                            %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Compute gradient smoothing         
                obj.avgGrads(:,pind) = obj.avgConst .* obj.grads(:,pind) + (1 - obj.avgConst) .* obj.avgGrads(:,pind);
                
                % Keep track of average gradient norm
                obj.gradNorm(1,pind)    = norm(obj.grads(:,pind))^2;
                obj.avgGradNorm(1,pind) = obj.avgConst .* abs(obj.gradNorm(1,pind)) + (1 - obj.avgConst) .* obj.avgGradNorm(1,pind);
                
                % Update stepsize
                obj.stepSize(:,pind) = (obj.stepSizeInit(:,pind) .* (obj.noiseVar(1, pind) / obj.noiseVarMin)^2) ./ (1e-4 + obj.stepSizeConst .* (obj.avgGrads(:,pind).^2));
                
                % Check if stepsize exceeds limits
                 obj.stepSize = max(obj.stepSizeMin, min(obj.stepSizeMax, obj.stepSize));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                              %
                % --- Langevin Monte Carlo --- %
                %                              %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             

                % Adjust the temperature
                obj.temper = obj.tempConst * log10(1 + (tind));

                % Do Langevin updates on all params to get a new proposed point
                paramProp = obj.param(:,pind) - obj.stepSize(:,pind) .* obj.avgGrads(:,pind) + ...
                                                obj.bEnableLangevin * sqrt(2 * obj.stepSize(:,pind) ./ obj.temper) .* randn(obj.numParams,1);

                % obj.param(:,pind) = paramProp;
                funcRatio = exp(obj.temper * (currentFuncVal - obj.cpe{1,pind}.evalObjectiveFunc(paramProp)));
                alpha     = min(1, funcRatio);
                unifRand  = rand; % ~ Unif[0,1]
                if unifRand <= alpha
                    obj.param(:,pind) = paramProp;
                    obj.bAccept(1,pind) = true;
                else
                    obj.bAccept(1,pind) = false;
                end
                obj.saveProposedParams(:,pind, tind) = paramProp;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                %
                % --- Save Data for Analysis --- %
                %                                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
                % Fill analysis arrays
                obj.saveParams(:,pind, tind)    = obj.param(:,pind);
                obj.saveGrads(:, pind, tind)    = obj.grads(:,pind);
                obj.saveHessTrc(pind,tind)      = obj.absHessianTrace(1,pind);
                obj.saveStepSize(:, pind, tind) = obj.stepSize(:,pind); 
                obj.saveGradNorm(pind, tind)    = obj.gradNorm(1,pind);
                obj.saveAvgGradNorm(pind, tind) = obj.avgGradNorm(1,pind);
                obj.saveObjFunc(pind, tind)     = currentFuncVal;
                obj.savebAccept(pind, tind)     = obj.bAccept(1,pind);
                obj.saveNoiseVar(pind, tind)    = obj.noiseVar(1,pind);

            end % end numParticles

            %%%%%%%%%%%%%%%%%%%
            %                 %
            % --- Waitbar --- %
            %                 %
            %%%%%%%%%%%%%%%%%%%

            % Update waitbar and message
            if getappdata(wbar, 'canceling')
                disp('Simulation cancelled!')

                obj.bStopSim = true;
                delete(wbar);
                return
            end    
            if mod(tind, 10) == 0
                
                waitbar(tind / obj.numIterLmc, wbar, ['LMC Iteration: ', num2str(tind), '/', num2str(obj.numIterLmc)]);

                % Update Live Plot
                if obj.bDisplayPlots
                    updateLivePlots_SL(obj, tind, ax);
                end

                % Print current parameters
                disp('--- Current parameters ---');
                for pind = 1:obj.numParticles
                    disp(['------ Particle  ', num2str(pind), ' = ', num2str(obj.param(:,pind).'), ' ------']);
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
            [~, minObjLmcInd] = min(obj.saveObjFunc(pind, :));
            obj.param(:,pind) = obj.saveParams(:,pind, minObjLmcInd); % This is the starting point for the next round of LMC updates
        end

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

        fprintf('\n');
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(['The optimum chirp parameters are =  ', num2str(obj.optParam.')]);
        disp(['Found by particle ', num2str(obj.bestParticleInd)]);
        fprintf('\n');

        % Finally, compute the scalar gains for CPE
        obj.cpe{1,obj.bestParticleInd} = obj.cpe{1,obj.bestParticleInd}.compScalarGains(obj.optParam);   

        % Store best particle index
        % Compute b vector (amp and phase too)

        % Housekeeping
        disp('Simulation complete!')
        delete(wbar);

    catch me
        delete(wbar); % Close wait bar when simulation errors out
        throw(me);
    end

end