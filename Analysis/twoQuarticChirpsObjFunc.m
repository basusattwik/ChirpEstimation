% MATLAB Code for Signal Analysis
clc
% close all
clearvars

phi{1,1} = [0, 10, 40, -70, 110].'; 
phi{1,2} = [0, 50, 60, -90, 105].';  
% Parameters
b = 1.3;
A1_true  = 1.0;
f01_true = 30;
k1_true  = 70.0;
c1_true  = 0;%-70;
d1_true  = 0;

A2_true  = 1;
f02_true = 45;
k2_true  = 85;
c2_true  = 0;%-90;
d2_true  = 0;

Ngs = 100;
N   = 250;
fs  = 1000;
T   = 1;  % Extended time duration
t   = 0:1/fs:T-1/fs; %linspace(0, T, 1000);  % Time vector

lim = 6;

% SNR increase factor
SNR_increase_factor = 1;

% Generate the signal without damping
s_measured_no_damping = SNR_increase_factor * (...
    A1_true * exp(-b * t) .* exp(1i * 2 * pi * (f01_true * t + k1_true * t.^2 + c1_true * t.^3 + d1_true * t.^4)) + ...
    A2_true * exp(-b * t) .* exp(1i * 2 * pi * (f02_true * t + k2_true * t.^2 + c2_true * t.^3 + d2_true * t.^4)));

s_measured_no_damping = s_measured_no_damping; %.* hann(length(t), "symmetric").';

% Finer grid for parameter search
f01_range = linspace(f01_true - lim , f01_true + lim , N);
f02_range = linspace(f02_true - lim , f02_true + lim , N);
k1_range  = linspace(k1_true - lim , k1_true + lim , N);
k2_range  = linspace(k2_true - lim , k2_true + lim , N);
% c1_range  = linspace(c1_true - lim , c1_true + lim , N);
% c2_range  = linspace(c2_true - lim , c2_true + lim , N);
% d1_range  = linspace(d1_true - lim , d1_true + lim , N);
% d2_range  = linspace(d2_true - lim , d2_true + lim , N);

% Preallocate MSE matrices
mse_f = zeros(length(f01_range), length(f02_range));
mse_k = zeros(length(k1_range), length(k2_range));
% mse_c = zeros(length(c1_range), length(c2_range));
% mse_d = zeros(length(d1_range), length(d2_range));
mse_f_gs = zeros(length(f01_range), length(f02_range));
mse_k_gs = zeros(length(k1_range), length(k2_range));
% mse_c_gs = zeros(length(c1_range), length(c2_range));
% mse_d_gs = zeros(length(d1_range), length(d2_range));

% Compute MSE for varying f01 and f02
temp = 0;
for i = 1:length(f01_range)
    for j = 1:length(f02_range)
        s_model = SNR_increase_factor * (...
            A1_true * exp(1i * 2 * pi * (f01_range(i) * t + k1_true * t.^2 + c1_true * t.^3 + d1_true * t.^4)) + ...
            A2_true * exp(1i * 2 * pi * (f02_range(j) * t + k2_true * t.^2 + c2_true * t.^3 + d2_true * t.^4)));
        mse_f(i, j) = sum(abs(s_measured_no_damping - s_model).^2);

        % for n = 1:Ngs
        %     s_model_gs = SNR_increase_factor * (...
        %         A1_true * exp(1i * 2 * pi * ((f01_range(i) + sig_gs*randn(1,1)) * t + k1_true * t.^2 + c1_true * t.^3)) + ...
        %         A2_true * exp(1i * 2 * pi * ((f02_range(j) + sig_gs*randn(1,1)) * t + k2_true * t.^2 + c2_true * t.^3)));
        %     temp = temp + mean(abs(s_measured_no_damping - s_model_gs).^2);
        % end
        % mse_f_gs(i, j) = temp/Ngs;
        % temp = 0;
    end
end

% Compute MSE for varying k1 and k2
for i = 1:length(k1_range)
    for j = 1:length(k2_range)
        s_model = SNR_increase_factor * (...
            A1_true * exp(1i * 2 * pi * (f01_true * t + k1_range(i) * t.^2 + c1_true * t.^3 + d1_true * t.^4)) + ...
            A2_true * exp(1i * 2 * pi * (f02_true * t + k2_range(j) * t.^2 + c2_true * t.^3 + d2_true * t.^4)));
        mse_k(i, j) = sum(abs(s_measured_no_damping - s_model).^2);

        % for n = 1:Ngs
        %     s_model_gs = SNR_increase_factor * (...
        %         A1_true * exp(1i * 2 * pi * (f01_true * t + (k1_range(i) + sig_gs*randn(1,1)) * t.^2 + c1_true * t.^3 + d1_true * t.^4)) + ...
        %         A2_true * exp(1i * 2 * pi * (f02_true * t + (k2_range(j) + sig_gs*randn(1,1)) * t.^2 + c2_true * t.^3 + d2_true * t.^4)));
        %     temp = temp + mean(abs(s_measured_no_damping - s_model_gs).^2);
        % end
        % mse_k_gs(i, j) = temp/Ngs;
        % temp = 0;

    end
end

% % Compute MSE for varying c1 and c2
% for i = 1:length(c1_range)
%     for j = 1:length(c2_range)
%         s_model = SNR_increase_factor * (...
%             A1_true * exp(1i * 2 * pi * (f01_true * t + k1_true * t.^2 + c1_range(i) * t.^3 + d1_true * t.^4)) + ...
%             A2_true * exp(1i * 2 * pi * (f02_true * t + k2_true * t.^2 + c2_range(j) * t.^3 + d2_true * t.^4)));
%         mse_c(i, j) = sum(abs(s_measured_no_damping - s_model).^2);
% 
%         % for n = 1:Ngs
%         %     s_model_gs = SNR_increase_factor * (...
%         %         A1_true * exp(1i * 2 * pi * (f01_true * t + k1_true * t.^2 + (c1_range(i) + sig_gs*randn(1,1)) * t.^3 + d1_true * t.^4)) + ...
%         %         A2_true * exp(1i * 2 * pi * (f02_true * t + k2_true * t.^2 + (c2_range(j) + sig_gs*randn(1,1)) * t.^3 + d2_true * t.^4)));
%         %     temp = temp + mean(abs(s_measured_no_damping - s_model_gs).^2);
%         % end
%         % mse_c_gs(i, j) = temp/Ngs;
%         % temp = 0;
%     end
% end
% 
% % Compute MSE for varying c1 and c2
% for i = 1:length(d1_range)
%     for j = 1:length(d2_range)
%         s_model = SNR_increase_factor * (...
%             A1_true * exp(1i * 2 * pi * (f01_true * t + k1_true * t.^2 + c1_true * t.^3 + d1_range(i) * t.^4)) + ...
%             A2_true * exp(1i * 2 * pi * (f02_true * t + k2_true * t.^2 + c2_true * t.^3 + d2_range(j) * t.^4)));
%         mse_d(i, j) = sum(abs(s_measured_no_damping - s_model).^2);
% 
%         % for n = 1:Ngs
%         %     s_model_gs = SNR_increase_factor * (...
%         %         A1_true * exp(1i * 2 * pi * (f01_true * t + k1_true * t.^2 + c1_true * t.^3 + (d1__range(i) + sig_gs*randn(1,1) * t.^4)) + ...
%         %         A2_true * exp(1i * 2 * pi * (f02_true * t + k2_true * t.^2 + c2_true + sig_gs*randn(1,1)) * t.^3 + (d2__range(j) + sig_gs*randn(1,1)) * t.^4)));
%         %     temp = temp + mean(abs(s_measured_no_damping - s_model_gs).^2); 
%         % end
%         % mse_c_gs(i, j) = temp/Ngs;
%         % temp = 0;
%     end
% end

%% Save data

% saveObjFuncData.freq = mse.
%% Plots
% 2D Plots for MSE

% Plot MSE for f01 and f02
% figure('windowstyle','docked');
% tiledlayout flow
% 
% nexttile
% contourf(f01_range, f02_range, mse_f', 10);
% colorbar;
% hold on;
% plot(f01_true, f02_true, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('Frequency f_{01}');
% ylabel('Frequency f_{02}');
% title('MSE for Varying f_{01} and f_{02} (No Damping)');
% 
% % Plot MSE for k1 and k2
% nexttile
% contourf(k1_range, k2_range, mse_k', 10);
% colorbar;
% hold on;
% plot(k1_true, k2_true, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('Chirp Rate k_1');
% ylabel('Chirp Rate k_2');
% title('MSE for Varying k_1 and k_2 (No Damping)');
% 
% % Plot MSE for c1 and c2
% nexttile
% contourf(c1_range, c2_range, mse_c', 10);
% colorbar;
% hold on;
% plot(c1_true, c2_true, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('Cubic Phase Coefficient c_1');
% ylabel('Cubic Phase Coefficient c_2');
% title('MSE for Varying c_1 and c_2 (No Damping)');
% 
% % Plot MSE for d1 and d2
% nexttile
% contourf(d1_range, d2_range, mse_d', 10);
% colorbar;
% hold on;
% plot(d1_true, d2_true, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('Quartic Phase Coefficient d_1');
% ylabel('Quartic Phase Coefficient d_2');
% title('MSE for Varying d_1 and d_2 (No Damping)');

% 3D Plots for MSE

% 3D Plot MSE for f01 and f02
fname = 'myFig';
hfig = figure;%('windowstyle','docked');
tiledlayout(1,1)

% nexttile
% surf(f01_range, f02_range, mse_f', 'EdgeColor', 'none');
% hold on;
% contour(f01_range, f02_range, mse_f', 10);
% scatter3(f01_true, f02_true, min(mse_f(:)), 100, 'r', 'filled');
% xlabel('$\varphi_{11}$', 'Interpreter','latex');
% ylabel('$\varphi_{21}$', 'Interpreter','latex');
zlabel('$J$', 'Interpreter','latex');
% zlim([min(mse_f(:))-5, max(mse_f(:)) + 5]);
% title('Ob Varying f_{01} and f_{02}');

% 3D Plot MSE for k1 and k2
nexttile
surf(k1_range, k2_range, mse_k', 'EdgeColor', 'none');
hold on;
contour(k1_range, k2_range, mse_k', 10);
scatter3(k1_true, k2_true, min(mse_k(:)), 100, 'r', 'filled');
xlabel('$\varphi_{12}$', 'Interpreter','latex');
ylabel('$\varphi_{22}$', 'Interpreter','latex');
zlabel('$J$', 'Interpreter','latex');
% zlim([min(mse_k(:))-5, max(mse_k(:)) + 5]);
% title('3D MSE for Varying k_1 and k_2 (No Damping)');

% % 3D Plot MSE for c1 and c2
% subplot(1,3,3); %nexttile
% surf(c1_range, c2_range, mse_c', 'EdgeColor', 'none');
% hold on;
% contour(c1_range, c2_range, mse_c', 10);
% scatter3(c1_true, c2_true, min(mse_c(:)), 100, 'r', 'filled');
% % colorbar
% xlabel('$\varphi_{13}$', 'Interpreter','latex');
% ylabel('$\varphi_{23}$', 'Interpreter','latex');
% zlabel('J');
% % title('3D MSE for Varying c_1 and c_2 (No Damping)');

set(findall(hfig,'-property','FontSize'),'FontSize',20,'FontName', 'Times')

% 3D Plot MSE for d1 and d2
% nexttile
% surf(d1_range, d2_range, mse_d', 'EdgeColor', 'none');
% hold on;
% contour(d1_range, d2_range, mse_d', 10);
% scatter3(d1_true, d2_true, min(mse_d(:)), 100, 'r', 'filled');
% xlabel('$\varphi_{14}$');
% ylabel('$\varphi_{24}$');
% zlabel('J');
% title('3D MSE for Varying d_1 and d_2 (No Damping)');

%3D Plots for Gibbs

% beta = 1;
% % 3D Plot MSE for f01 and f02
% figure('windowstyle','docked');
% tiledlayout flow
% 
% nexttile
% surf(f01_range, f02_range, exp(-beta * mse_f'), 'EdgeColor', 'none');
% hold on;
% scatter3(f01_true, f02_true, min(mse_f(:)), 100, 'r', 'filled');
% xlabel('Frequency f_{01}');
% ylabel('Frequency f_{02}');
% zlabel('J');
% title('3D MSE for Varying f_{01} and f_{02} (No Damping)');
% 
% % 3D Plot MSE for k1 and k2
% nexttile
% surf(k1_range, k2_range, exp(-beta * mse_k'), 'EdgeColor', 'none');
% hold on;
% scatter3(k1_true, k2_true, min(mse_k(:)), 100, 'r', 'filled');
% xlabel('Chirp Rate k_1');
% ylabel('Chirp Rate k_2');
% zlabel('J');
% title('3D MSE for Varying k_1 and k_2 (No Damping)');
% 
% % 3D Plot Gibbs for c1 and c2
% nexttile
% surf(c1_range, c2_range, exp(-beta * mse_c'), 'EdgeColor', 'none');
% hold on;
% scatter3(c1_true, c2_true, min(mse_c(:)), 100, 'r', 'filled');
% xlabel('Cubic Phase Coefficient c_1');
% ylabel('Cubic Phase Coefficient c_2');
% zlabel('J');
% title('3D MSE for Varying c_1 and c_2 (No Damping)');
% 
% % 3D Plot Gibbs for c1 and c2
% nexttile
% surf(d1_range, d2_range, exp(-beta * mse_d'), 'EdgeColor', 'none');
% hold on;
% scatter3(d1_true, d2_true, min(mse_d(:)), 100, 'r', 'filled');
% xlabel('Quartic Phase Coefficient d_1');
% ylabel('Quartic Phase Coefficient d_2');
% zlabel('J');
% title('3D MSE for Varying d_1 and d_2 (No Damping)');


% makeNiceFigure(hfig, fname)
%%

function [] = makeNiceFigure(hfig, fname)

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.2; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize', 20, 'FontName', 'Times') % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig, fname,'-dpng','-painters')

end