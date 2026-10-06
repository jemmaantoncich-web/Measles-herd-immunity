% BIOMATHS 374 - MEASLES HERD IMMUNITY PROJECT
% Two-patch SEIR-V metapopulation model: final results

%
% The model equations are identical to the midyear version (equations 1-8).
% Requires MATLAB R2020a or later (uses yline, table, bar CData and exportgraphics).

function measles_model()

% =====================================================================
% 1. BASE PARAMETERS (Table 3)
% =====================================================================
 p.beta  = 1.25;       % Transmission rate (day^-1)
 p.sigma = 0.125;      % Progression rate E -> I (day^-1), 8-day incubation
 p.gamma = 0.125;      % Recovery rate I -> R (day^-1), 8-day infectious period
 p.mu    = 0.0000418;  % Natural birth/death rate (day^-1), 1/(65*365)
 p.eps1  = 0.85;       % Efficacy of dose 1
 p.eps2  = 0.97;       % Efficacy of dose 2
 p.N1    = 400000;     % Population of Patch 1 (illustrative)
 p.N2    = 1600000;    % Population of Patch 2 (illustrative)
 p.phi   = 0.1;        % Mixing fraction
 tEnd    = 730;        % Length of each simulation (days)


 col.blue   = [0.00 0.45 0.74];   % susceptible / scenario 2
 col.green  = [0.47 0.67 0.19];   % vaccinated
 col.yellow = [0.93 0.69 0.13];   % exposed
 col.red    = [0.85 0.33 0.10];   % infectious
 col.purple = [0.49 0.18 0.56];   % recovered / scenario 3

 R0 = calc_R0(p);
 fprintf('Basic reproduction number R0 = %.4f\n', R0);

% Transmission rate that gives R0 = 20 (used for the R0 = 20 check)
 beta20 = 20 * (p.sigma + p.mu) * (p.gamma + p.mu) / p.sigma;
 fprintf('Transmission rate for R0 = 20: beta = %.4f\n\n', beta20);

% =====================================================================
% 2. SCENARIOS (Table 4)
% =====================================================================
% Scenario 1: uniform 95% coverage in both patches
 p1 = p;
 p1.v1_1 = 0.05;   p1.v2_1 = 0.90;    % Patch 1: 95% total
 p1.v1_2 = 0.05;   p1.v2_2 = 0.90;    % Patch 2: 95% total

% Scenario 2: clustered coverage, same 95% population-weighted average
 p2 = p;
 p2.v1_1 = 0.20;   p2.v2_1 = 0.60;    % Patch 1: 80% total
 p2.v1_2 = 0.0475; p2.v2_2 = 0.94;    % Patch 2: 98.75% total

% Scenario 3: Capricorn (Patch 1) and City of Johannesburg (Patch 2)
% v2 = second-dose (MCV2) coverage; v1 = MCV1 - MCV2 (one dose only)
 p3 = p;
 p3.N1   = 1447103;  p3.N2 = 4803262;  % Census 2022 populations
 p3.v1_1 = 0.081;  p3.v2_1 = 0.755;    % Capricorn: MCV1 83.6%, MCV2 75.5%
 p3.v1_2 = 0.045;  p3.v2_2 = 0.904;    % Johannesburg: MCV1 94.9%, MCV2 90.4%

 scen  = {p1, p2, p3};
 names = {'Scenario 1: Uniform', 'Scenario 2: Clustered', ...
          'Scenario 3: Capricorn-Johannesburg'};
 letters = {'(a)', '(b)', '(c)', '(d)', '(e)'};

% =====================================================================
% 3. RUN ALL SCENARIOS (R0 = 10 and R0 = 20)
% =====================================================================
 nS = numel(scen);
 Rc = zeros(nS,1); peakI = zeros(nS,1); peakDay = zeros(nS,1); endDay = zeros(nS,1);
 totInf = zeros(nS,1); totInf1 = zeros(nS,1); totInf2 = zeros(nS,1);
 Rc20 = zeros(nS,1); peakI20 = zeros(nS,1); totInf20 = zeros(nS,1);
 T = cell(nS,1); Y = cell(nS,1); T20 = cell(nS,1); Y20 = cell(nS,1);

 for k = 1:nS
   % Main results, R0 = 10
   [T{k}, Y{k}, Rc(k)] = run_simulation(scen{k}, tEnd);
   [peakI(k), peakDay(k), endDay(k), totInf(k), totInf1(k), totInf2(k)] = ...
       outbreak_size(T{k}, Y{k}, scen{k});

   % Check with R0 = 20 (only beta changes)
   q = scen{k};  q.beta = beta20;
   [T20{k}, Y20{k}, Rc20(k)] = run_simulation(q, tEnd);
   [peakI20(k), ~, ~, totInf20(k)] = outbreak_size(T20{k}, Y20{k}, q);
 end

% Summary table (Table 5 in the report)
 Scenario = names';
 summary = table(Scenario, Rc, peakI, peakDay, endDay, totInf1, totInf2, totInf, ...
   'VariableNames', {'Scenario','Rc','PeakInfectious','PeakDay','OutbreakEndDay', ...
                     'TotalInfected_Patch1','TotalInfected_Patch2','TotalInfected'});
 disp('================== TABLE 5: SUMMARY OF SCENARIOS (R0 = 10) ==================');
 disp(summary);
 disp('(OutbreakEndDay = first day after the peak with I1 + I2 < 1; NaN = no outbreak or not ended)');
 fprintf('\n');

 summary20 = table(Scenario, Rc, Rc20, peakI20, totInf20, Rc20 ./ Rc, ...
   'VariableNames', {'Scenario','Rc_R0_10','Rc_R0_20','PeakInfectious_R0_20', ...
                     'TotalInfected_R0_20','Ratio_Rc'});
 disp('================== R0 = 20 CHECK ==================');
 disp(summary20);
 fprintf('Clustered/uniform Rc ratio: %.2f (R0 = 10) and %.2f (R0 = 20)\n\n', ...
         Rc(2)/Rc(1), Rc20(2)/Rc20(1));

% =====================================================================
% 4. FIGURE 2: INFECTIOUS PROPORTION, ALL SCENARIOS (main text)
% =====================================================================
% Each panel has its own y-axis because Scenario 1 values are tiny.
% MATLAB shows the scale (e.g. x10^-6) above the y-axis; mention it in the
% results paragraph under the figure.
 plot_infectious(T, Y, scen, names, Rc, letters, col.red, ...
                 'Figure 2: Infectious proportion', 'Figure2.png');
 fprintf('Largest infectious proportion in Scenario 1: Patch 1 %.2e, Patch 2 %.2e\n\n', ...
         max(Y{1}(:,5))/p1.N1, max(Y{1}(:,11))/p1.N2);

% =====================================================================
% 4b. FIGURE 3 AND TABLE 6: SUSCEPTIBLE AND VACCINATED (main text)
% =====================================================================
% Figure 3: (a) S and (b) V1 + V2 for Scenario 2, (c) S and (d) V1 + V2
% for Scenario 3. Scenario 1 is left out because S and V change only in
% the sixth decimal place (see Figure A1).
 fig3 = figure('Name','Figure 3: Susceptible and vaccinated','Color','w', ...
               'Position',[150 150 900 700]);
 panel = 0;
 for k = 2:3
   t = T{k}; Yk = Y{k}; pk = scen{k};
   Sp = {Yk(:,1)/pk.N1, Yk(:,7)/pk.N2};                       % S1, S2
   Vp = {(Yk(:,2) + Yk(:,3))/pk.N1, (Yk(:,8) + Yk(:,9))/pk.N2}; % V1 + V2
   data = {Sp, Vp};
   lbl  = {'Susceptible (S)', 'Vaccinated (V_1 + V_2)'};
   cl   = {col.blue, col.green};
   for c = 1:2
     panel = panel + 1;
     subplot(2, 2, panel);
     plot(t, data{c}{1}, '-',  'Color', cl{c}, 'LineWidth', 1.8); hold on;
     plot(t, data{c}{2}, '--', 'Color', cl{c}, 'LineWidth', 1.8); hold off;
     grid on; xlim([0 t(end)]);
     xlabel('Time (days)'); ylabel('Proportion of patch');
     title({sprintf('%s %s', letters{panel}, names{k}), lbl{c}});
     legend('Patch 1', 'Patch 2', 'Location', 'best');
   end
 end
 exportgraphics(fig3, 'Figure3.png', 'Resolution', 300);

% Values to quote under Figure 3
 for k = 2:3
   pk = scen{k}; Yk = Y{k}; t = T{k};
   for j = 1:2
     o = 6*(j - 1); Nj = pk.(sprintf('N%d', j));
     [Smin, iS] = min(Yk(:,o+1));
     fprintf('%s, Patch %d: S %.3f -> min %.3f (day %.1f) -> %.3f at day %d; V %.3f -> %.3f\n', ...
             names{k}, j, Yk(1,o+1)/Nj, Smin/Nj, t(iS), Yk(end,o+1)/Nj, t(end), ...
             (Yk(1,o+2) + Yk(1,o+3))/Nj, (Yk(end,o+2) + Yk(end,o+3))/Nj);
   end
 end
 fprintf('\n');

% Table 6: where the infections in each patch came from
 rowNames = {}; fromS = []; fromV1 = []; fromV2 = [];
 for k = 2:3
   B = infections_by_status(T{k}, Y{k}, scen{k});   % 2 x 3: patches x (S, V1, V2)
   for j = 1:2
     rowNames{end+1,1} = sprintf('%s, Patch %d', names{k}, j); %#ok<AGROW>
     fromS(end+1,1) = B(j,1); fromV1(end+1,1) = B(j,2); fromV2(end+1,1) = B(j,3); %#ok<AGROW>
   end
 end
 tot = fromS + fromV1 + fromV2;
 tab6 = table(rowNames, round(fromS), round(fromV1), round(fromV2), round(tot), ...
              100*fromS./tot, 100*fromV1./tot, 100*fromV2./tot, ...
   'VariableNames', {'Patch','From_S','From_V1','From_V2','Total', ...
                     'Pct_S','Pct_V1','Pct_V2'});
 disp('================== TABLE 6: INFECTIONS BY VACCINATION STATUS ==================');
 disp(tab6);
 disp('(Totals exclude the one initially exposed person in Patch 1)');
 fprintf('\n');

% =====================================================================
% 5. FIGURE 4: Rc AGAINST THE MIXING FRACTION (main text)
% =====================================================================
% phi is a fraction of contacts, so it must lie between 0 and 1.
 phi_values = linspace(0, 1, 101);
 Rc_phi  = zeros(size(phi_values));   % Scenario 2: clustered
 Rc_phi3 = zeros(size(phi_values));   % Scenario 3: Capricorn-Johannesburg
 for i = 1:numel(phi_values)
   pt = p2;  pt.phi = phi_values(i);  Rc_phi(i)  = calc_Rc(pt);
   pt = p3;  pt.phi = phi_values(i);  Rc_phi3(i) = calc_Rc(pt);
 end
 fprintf('Clustered: Rc = %.4f at phi = 0, %.4f at phi = 0.5, %.4f at phi = 1\n', ...
         Rc_phi(1), Rc_phi(51), Rc_phi(end));
 fprintf('Capricorn-Johannesburg: Rc = %.4f at phi = 0, %.4f at phi = 0.5, %.4f at phi = 1\n', ...
         Rc_phi3(1), Rc_phi3(51), Rc_phi3(end));
 fprintf('Minimum Rc over 0 <= phi <= 1: clustered %.4f, district pair %.4f\n\n', ...
         min(Rc_phi), min(Rc_phi3));

 fig4 = figure('Name','Figure 4: Rc vs phi','Color','w','Position',[200 200 600 400]);
 plot(phi_values, Rc_phi,  '-', 'Color', col.blue,   'LineWidth', 2); hold on;
 plot(phi_values, Rc_phi3, '-', 'Color', col.purple, 'LineWidth', 2);
 yline(1, '--', 'R_c = 1', 'Color', 'k', 'LineWidth', 1.5, ...
       'LabelHorizontalAlignment', 'left');
 hold off; grid on; xlim([0 1]); ylim([0.8 2.6]);
 xlabel('Cross-patch mixing fraction \phi (dimensionless)');
 ylabel('Control reproduction number R_c');
 legend('Scenario 2: Clustered', 'Scenario 3: Capricorn-Johannesburg', 'Location', 'northeast');
 exportgraphics(fig4, 'Figure4.png', 'Resolution', 300);

% =====================================================================
% 6. FIGURE 5: SENSITIVITY ANALYSIS, SCENARIOS 2 AND 3 (main text)
% =====================================================================
% Normalised forward sensitivity index Upsilon_p = (dRc/dp)*(p/Rc)
% (equation 13), calculated with a central difference (+/- 1% change in p).
 pnames  = {'beta','sigma','gamma','eps1','eps2','phi','v1_1','v2_1'};
 plabels = {'\beta','\sigma','\gamma','\epsilon_1','\epsilon_2','\phi','v_{1,1}','v_{2,1}'};
 Ups2 = sensitivity(p2, pnames);
 Ups3 = sensitivity(p3, pnames);

 Parameter = pnames';
 disp('============ SENSITIVITY INDICES ============');
 disp(table(Parameter, Ups2, Ups3, 'VariableNames', ...
      {'Parameter','Scenario2_Clustered','Scenario3_District'}));

 fig5 = figure('Name','Figure 5: Sensitivity indices','Color','w','Position',[250 250 1000 400]);
 yl = [min([Ups2; Ups3]) - 0.3, max([Ups2; Ups3]) + 0.3];   % same y-axis for both panels
 UpsAll = {Ups2, Ups3};
 sNames = {names{2}, names{3}};
 for s = 1:2
   subplot(1, 2, s);
   Ups = UpsAll{s};
   b = bar(Ups, 'FaceColor', 'flat', 'EdgeColor', 'k');
   for i = 1:numel(Ups)              % red = increases Rc, blue = decreases Rc
     if Ups(i) >= 0
       b.CData(i,:) = col.red;
     else
       b.CData(i,:) = col.blue;
     end
   end
   set(gca, 'XTick', 1:numel(pnames), 'XTickLabel', plabels);
   ylim(yl); grid on;
   ylabel('Sensitivity index \Upsilon_p');
   title(sprintf('%s %s', letters{s}, sNames{s}));
 end
 exportgraphics(fig5, 'Figure5.png', 'Resolution', 300);

% =====================================================================
% 7. CRITICAL COVERAGE (Scenario 3)
% =====================================================================
% Smallest second-dose coverage v2_1 in Capricorn (Patch 1) that gives
% Rc = 1, keeping one-dose-only coverage v1_1 and Patch 2 fixed.
 f = @(v) calc_Rc(setfield(p3, 'v2_1', v)) - 1; %#ok<SFLD>
 vmin = p3.v2_1;
 vmax = 1 - p3.v1_1;                % total coverage cannot exceed 1
 fprintf('\n============ CRITICAL COVERAGE (Capricorn) ============\n');
 if f(vmin) <= 0
   fprintf('Rc is already below 1 at current coverage.\n\n');
 elseif f(vmax) > 0
   fprintf('Rc stays above 1 even at maximum coverage; no critical value exists.\n\n');
 else
   v_crit = fzero(f, [vmin vmax]);
   fprintf('Current MCV2 coverage: %.1f%%\n', 100*p3.v2_1);
   fprintf('Critical MCV2 coverage for Rc = 1: %.1f%% (increase of %.1f percentage points)\n\n', ...
           100*v_crit, 100*(v_crit - p3.v2_1));
 end

% =====================================================================
% 8. APPENDIX FIGURES
% =====================================================================
% Figures A1-A3: every compartment for one scenario, panels (a)-(e).
 cmpNames = {'Susceptible (S)', 'Vaccinated (V_1 + V_2)', 'Exposed (E)', ...
             'Infectious (I)', 'Recovered (R)'};
 cmpCols  = {col.blue, col.green, col.yellow, col.red, col.purple};
 for k = 1:nS
   plot_all_compartments(T{k}, Y{k}, scen{k}, names{k}, Rc(k), cmpNames, cmpCols, ...
                         letters, sprintf('Figure A%d: %s', k, names{k}), ...
                         sprintf('FigureA%d.png', k));
 end

% Figure A4: infectious proportion for all scenarios with R0 = 20
 plot_infectious(T20, Y20, scen, names, Rc20, letters, col.red, ...
                 'Figure A4: Infectious proportion, R0 = 20', 'FigureA4.png');

 fprintf('Main figures saved as Figure2.png to Figure5.png\n');
 fprintf('Appendix figures saved as FigureA1.png to FigureA4.png\n');
end

% -------------------------------------------------------------------------
% Basic reproduction number (no vaccination)
% -------------------------------------------------------------------------
function R0 = calc_R0(p)
 R0 = p.beta*p.sigma / ((p.sigma + p.mu)*(p.gamma + p.mu));
end

% -------------------------------------------------------------------------
% Next generation matrix and Rc (equations 11 and 12)
% -------------------------------------------------------------------------
function R_c = calc_Rc(p)
% Effective susceptibles at the disease-free equilibrium
 S_eff_1 = p.N1 * ((1 - p.v1_1 - p.v2_1) + (1 - p.eps1)*p.v1_1 + (1 - p.eps2)*p.v2_1);
 S_eff_2 = p.N2 * ((1 - p.v1_2 - p.v2_2) + (1 - p.eps1)*p.v1_2 + (1 - p.eps2)*p.v2_2);

% F: new infections (rows/columns ordered E1, I1, E2, I2)
 F = [0, p.beta*(1 - p.phi)*(S_eff_1/p.N1), 0, p.beta*p.phi*(S_eff_1/p.N2);
      0, 0, 0, 0;
      0, p.beta*p.phi*(S_eff_2/p.N1), 0, p.beta*(1 - p.phi)*(S_eff_2/p.N2);
      0, 0, 0, 0];

% V: transitions out of the infected compartments
 V = [p.sigma + p.mu, 0, 0, 0;
      -p.sigma, p.gamma + p.mu, 0, 0;
      0, 0, p.sigma + p.mu, 0;
      0, 0, -p.sigma, p.gamma + p.mu];

 K = F / V;                          % next generation matrix K = F V^-1
 R_c = max(abs(eig(K)));             % spectral radius
end

% -------------------------------------------------------------------------
% Normalised sensitivity indices of Rc (equation 13)
% -------------------------------------------------------------------------
function Ups = sensitivity(p, pnames)
 h = 0.01;
 Rc_base = calc_Rc(p);
 Ups = zeros(numel(pnames),1);
 for i = 1:numel(pnames)
   pu = p; pd = p;
   pu.(pnames{i}) = p.(pnames{i}) * (1 + h);
   pd.(pnames{i}) = p.(pnames{i}) * (1 - h);
   dRc = (calc_Rc(pu) - calc_Rc(pd)) / (2*h*p.(pnames{i}));
   Ups(i) = dRc * p.(pnames{i}) / Rc_base;
 end
end

% -------------------------------------------------------------------------
% Solve the ODE system for one scenario
% -------------------------------------------------------------------------
function [t, Y, R_c] = run_simulation(p, tEnd)
 R_c = calc_Rc(p);

% Initial conditions: one exposed individual in Patch 1, everyone else at
% the disease-free equilibrium proportions (equation 9)
 E1_0 = 1; I1_0 = 0; R1_0 = 0;
 S1_0   = (p.N1 - E1_0) * (1 - p.v1_1 - p.v2_1);
 V1_1_0 = (p.N1 - E1_0) * p.v1_1;
 V2_1_0 = (p.N1 - E1_0) * p.v2_1;
 E2_0 = 0; I2_0 = 0; R2_0 = 0;
 S2_0   = p.N2 * (1 - p.v1_2 - p.v2_2);
 V1_2_0 = p.N2 * p.v1_2;
 V2_2_0 = p.N2 * p.v2_2;
 Y0 = [S1_0; V1_1_0; V2_1_0; E1_0; I1_0; R1_0; S2_0; V1_2_0; V2_2_0; E2_0; I2_0; R2_0];

 tspan = 0:0.5:tEnd;                 % output every half day for smooth plots
 options = odeset('NonNegative', 1:12, 'RelTol', 1e-8, 'AbsTol', 1e-6);
 [t, Y] = ode45(@(t,Y) seir_v_odes(t, Y, p), tspan, Y0, options);
end

% -------------------------------------------------------------------------
% Outbreak size measures (Section 3.10)
% -------------------------------------------------------------------------
function [peakI, peakDay, endDay, totInf, totInf1, totInf2] = outbreak_size(t, Y, p)
% Peak number of people infectious at the same time (I1 + I2)
 Itot = Y(:,5) + Y(:,11);
 [peakI, idx] = max(Itot);
 peakDay = t(idx);

% Outbreak end: first time after the peak that fewer than one person is
% infectious. NaN if there was no outbreak (peak below 1) or it has not ended.
 after = find(t > peakDay & Itot < 1, 1);
 if peakI < 1 || isempty(after)
   endDay = NaN;
 else
   endDay = t(after);
 end

% Total infected = cumulative incidence: all new infections over the run,
% i.e. the integral of lambda_j*(S_j + (1-eps1)V1_j + (1-eps2)V2_j) dt,
% plus the one initially exposed person. This only post-processes the
% solution; the model equations are unchanged.
 I1 = Y(:,5); I2 = Y(:,11);
 lambda_1 = p.beta * ((1 - p.phi) * I1 / p.N1 + p.phi * I2 / p.N2);
 lambda_2 = p.beta * ((1 - p.phi) * I2 / p.N2 + p.phi * I1 / p.N1);
 newInf1 = lambda_1 .* (Y(:,1) + (1 - p.eps1)*Y(:,2) + (1 - p.eps2)*Y(:,3));
 newInf2 = lambda_2 .* (Y(:,7) + (1 - p.eps1)*Y(:,8) + (1 - p.eps2)*Y(:,9));
 totInf1 = trapz(t, newInf1) + Y(1,4);   % + initial exposed person
 totInf2 = trapz(t, newInf2);
 totInf  = totInf1 + totInf2;
end

% -------------------------------------------------------------------------
% Cumulative infections in each patch split by vaccination status.
% Row j = patch j; columns = infections of unvaccinated (S), one-dose (V1)
% and two-dose (V2) individuals: integral of lambda_j*S_j,
% lambda_j*(1-eps1)*V1_j and lambda_j*(1-eps2)*V2_j over the simulation.
% -------------------------------------------------------------------------
function B = infections_by_status(t, Y, p)
 I1 = Y(:,5); I2 = Y(:,11);
 lambda = {p.beta * ((1 - p.phi) * I1 / p.N1 + p.phi * I2 / p.N2), ...
           p.beta * ((1 - p.phi) * I2 / p.N2 + p.phi * I1 / p.N1)};
 B = zeros(2,3);
 for j = 1:2
   o = 6*(j - 1);
   B(j,1) = trapz(t, lambda{j} .* Y(:,o+1));
   B(j,2) = trapz(t, lambda{j} .* (1 - p.eps1) .* Y(:,o+2));
   B(j,3) = trapz(t, lambda{j} .* (1 - p.eps2) .* Y(:,o+3));
 end
end

% -------------------------------------------------------------------------
% Infectious proportion, one panel per scenario: (a), (b), (c).
% Patch 1 = solid line, Patch 2 = dashed line.
% -------------------------------------------------------------------------
function plot_infectious(Tc, Yc, Pc, Nc, Rcc, letters, lineCol, figName, fileName)
 nS = numel(Tc);
 fig = figure('Name', figName, 'Color', 'w', 'Position', [100 100 400*nS 380]);
 for k = 1:nS
   t = Tc{k}; Y = Yc{k}; p = Pc{k};
   subplot(1, nS, k);
   plot(t, Y(:,5)/p.N1,  '-',  'Color', lineCol, 'LineWidth', 1.8); hold on;
   plot(t, Y(:,11)/p.N2, '--', 'Color', lineCol, 'LineWidth', 1.8); hold off;
   grid on; xlim([0 t(end)]);
   xlabel('Time (days)');
   ylabel('Proportion of patch infectious');
   title({sprintf('%s %s', letters{k}, Nc{k}), sprintf('R_c = %.2f', Rcc(k))});
   legend('Patch 1', 'Patch 2', 'Location', 'best');
 end
 exportgraphics(fig, fileName, 'Resolution', 300);
end

% -------------------------------------------------------------------------
% All five compartments for one scenario: panels (a)-(e) (appendix).
% Patch 1 = solid line, Patch 2 = dashed line.
% -------------------------------------------------------------------------
function plot_all_compartments(t, Y, p, scenName, Rc, cmpNames, cmpCols, letters, figName, fileName)
 P1 = [Y(:,1), Y(:,2) + Y(:,3), Y(:,4), Y(:,5), Y(:,6)] / p.N1;
 P2 = [Y(:,7), Y(:,8) + Y(:,9), Y(:,10), Y(:,11), Y(:,12)] / p.N2;
 fig = figure('Name', figName, 'Color', 'w', 'Position', [100 100 1200 700]);
 for c = 1:5
   subplot(2, 3, c);
   plot(t, P1(:,c), '-',  'Color', cmpCols{c}, 'LineWidth', 1.8); hold on;
   plot(t, P2(:,c), '--', 'Color', cmpCols{c}, 'LineWidth', 1.8); hold off;
   grid on; xlim([0 t(end)]);
   xlabel('Time (days)');
   ylabel('Proportion of patch');
   title(sprintf('%s %s', letters{c}, cmpNames{c}));
   legend('Patch 1', 'Patch 2', 'Location', 'best');
 end
 % Empty sixth panel used for the scenario label
 subplot(2, 3, 6); axis off;
 text(0.05, 0.6, scenName, 'FontSize', 12, 'FontWeight', 'bold');
 text(0.05, 0.4, sprintf('R_c = %.2f', Rc), 'FontSize', 12);
 exportgraphics(fig, fileName, 'Resolution', 300);
end

% -------------------------------------------------------------------------
% ODE system (equations 1-8)
% -------------------------------------------------------------------------
function dYdt = seir_v_odes(~, Y, p)
 S1 = Y(1); V1_1 = Y(2); V2_1 = Y(3); E1 = Y(4); I1 = Y(5); R1 = Y(6);
 S2 = Y(7); V1_2 = Y(8); V2_2 = Y(9); E2 = Y(10); I2 = Y(11); R2 = Y(12);

% Force of infection (equations 1 and 2)
 lambda_1 = p.beta * ((1 - p.phi) * (I1 / p.N1) + p.phi * (I2 / p.N2));
 lambda_2 = p.beta * ((1 - p.phi) * (I2 / p.N2) + p.phi * (I1 / p.N1));

% Patch 1 (equations 3-8, j = 1)
 dS1_dt   = p.mu * p.N1 * (1 - p.v1_1 - p.v2_1) - lambda_1 * S1 - p.mu * S1;
 dV1_1_dt = p.mu * p.N1 * p.v1_1 - lambda_1 * (1 - p.eps1) * V1_1 - p.mu * V1_1;
 dV2_1_dt = p.mu * p.N1 * p.v2_1 - lambda_1 * (1 - p.eps2) * V2_1 - p.mu * V2_1;
 new_inf_1 = lambda_1 * (S1 + (1 - p.eps1) * V1_1 + (1 - p.eps2) * V2_1);
 dE1_dt = new_inf_1 - (p.sigma + p.mu) * E1;
 dI1_dt = p.sigma * E1 - (p.gamma + p.mu) * I1;
 dR1_dt = p.gamma * I1 - p.mu * R1;

% Patch 2 (equations 3-8, j = 2)
 dS2_dt   = p.mu * p.N2 * (1 - p.v1_2 - p.v2_2) - lambda_2 * S2 - p.mu * S2;
 dV1_2_dt = p.mu * p.N2 * p.v1_2 - lambda_2 * (1 - p.eps1) * V1_2 - p.mu * V1_2;
 dV2_2_dt = p.mu * p.N2 * p.v2_2 - lambda_2 * (1 - p.eps2) * V2_2 - p.mu * V2_2;
 new_inf_2 = lambda_2 * (S2 + (1 - p.eps1) * V1_2 + (1 - p.eps2) * V2_2);
 dE2_dt = new_inf_2 - (p.sigma + p.mu) * E2;
 dI2_dt = p.sigma * E2 - (p.gamma + p.mu) * I2;
 dR2_dt = p.gamma * I2 - p.mu * R2;

 dYdt = [dS1_dt; dV1_1_dt; dV2_1_dt; dE1_dt; dI1_dt; dR1_dt; ...
         dS2_dt; dV1_2_dt; dV2_2_dt; dE2_dt; dI2_dt; dR2_dt];
end
