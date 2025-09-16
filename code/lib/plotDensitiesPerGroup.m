close all
% Load data
T = readtable(fullfile('F:','lab','MouseMusclePOGLUT1',...
'results','numberOfPax7PerImage_17-Jul-2025.xlsx'),'VariableNamingRule','preserve');

%% Extract genotype and age from Row names
rowNames = T.Row; % adjust if your row names column has a different name
genotype = strings(size(rowNames));
age = strings(size(rowNames));

for i = 1:numel(rowNames)
    str = string(rowNames{i});
    if contains(str, "WT", 'IgnoreCase', true)
        genotype(i) = "WT";
    elseif contains(str, "KIKO", 'IgnoreCase', true)
        genotype(i) = "KIKO";
    else
        genotype(i) = "Unknown";
    end

    tokens = regexp(str, "p\d+", 'match');
    if ~isempty(tokens)
        age(i) = tokens{1};
    else
        age(i) = "NA";
    end
end

% Add to table
T.Genotype = categorical(genotype);
T.Age = categorical(age);

%% Define variables and plotting parameters
densityVars = ["whole muscle nuclei/um2", ...
               "nuclei/um2 g1","nuclei/um2 g2","nuclei/um2 g3", ...
               "nuclei/um2 g4","nuclei/um2 g5", ...
               "nuclei/um2 g7","nuclei/um2 g8","nuclei/um2 g9", ...
               "nuclei/um2 g10","nuclei/um2 g11","nuclei/um2 g12", ...
               "nuclei/um2 g13","nuclei/um2 g14","nuclei/um2 g15"];

agesOrder = {'p4','p7','p10','p15','p20','p35'};
genotypesList = ["WT","KIKO"];
colors = [0 0.4470 0.7410;   % WT blue
          0.8500 0.3250 0.0980]; % KIKO orange

%% Folder to save plots
saveFolder = fullfile('F:','lab','MouseMusclePOGLUT1','results','plots');
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end


% Plot each group
%% Loop over density variables
for v = 1:numel(densityVars)
    varName = densityVars(v);
     % Create figure: smaller size [width height]
    f = figure('Color','w','Units','pixels','Position',[100 100 800 600]); 
    hold on;
    title(strrep(varName,'_',' '), 'Interpreter','none', 'FontSize',14);


    % Plot points and error bars
    for a = 1:numel(agesOrder)
        for g = 1:numel(genotypesList)
            sel = (T.Genotype == genotypesList(g)) & (T.Age == agesOrder{a});
            vals = T{sel, varName};
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end

            x = (a-1)*3 + g; % WT left, KIKO right
            jitter = (rand(size(vals))-0.5)*0.25;

            % Scatter points
            scatter(x + jitter, vals, 45, ...
                'MarkerFaceColor', colors(g,:), ...
                'MarkerEdgeColor','k', ...
                'MarkerFaceAlpha',0.5, ...
                'MarkerEdgeAlpha',0.5);

            % Mean ± SD
            mu = mean(vals);
            sd = std(vals);
            errorbar(x, mu, sd, 'Color', colors(g,:), ...
                'LineWidth',1.8, 'CapSize',10);
        end

        % --- Statistical significance ---
        valsWT   = T{(T.Genotype=="WT") & (T.Age==agesOrder{a}), varName};
        valsKIKO = T{(T.Genotype=="KIKO") & (T.Age==agesOrder{a}), varName};
        valsWT = valsWT(~isnan(valsWT));
        valsKIKO = valsKIKO(~isnan(valsKIKO));
        if isempty(valsWT) || isempty(valsKIKO) || length(valsWT)<4 || length(valsKIKO)<4, continue; end

        [~,p] = compareMeans(valsWT, valsKIKO);

        if p < 0.001, stars = '***';
        elseif p < 0.01, stars = '**';
        elseif p < 0.05, stars = '*';
        else stars = 'n.s.';
        end

        x1 = (a-1)*3 + 1; % WT
        x2 = (a-1)*3 + 2; % KIKO
        ymax = max([valsWT; valsKIKO]);
        y = ymax + 0.1*ymax;

        if p < 0.05
            plot([x1 x2],[y y],'k-','LineWidth',1.2);
            text(mean([x1 x2]), y + 0.05*ymax, stars, ...
            'HorizontalAlignment','center','FontSize',12);
        end
        
    end

    % --- Axes formatting ---
    xticks(1.5:3:3*numel(agesOrder));
    xticklabels(agesOrder);
    ylabel('Nuclei density (\mum^{-2})','FontSize',12);
    set(gca,'FontSize',12,'LineWidth',1.2,'Box','off');

    legend({'','','','','','','WT','','KIKO'},'Location','best','Box','off');
    hold off;

     % --- Save figure ---
    saveName = fullfile(saveFolder, strrep(varName,'nuclei/um2','nucleiDensity'));

    % Use exportgraphics if MATLAB supports it
    try
        exportgraphics(f, strjoin([saveName,'.tif']), 'Resolution',300);
        exportgraphics(f, strjoin([saveName,'.png']), 'Resolution',300);
    catch
        % fallback to print for older MATLAB
        print(f, saveName, '-dtiff', '-r300');
        print(f, saveName, '-dpng', '-r300');
    end

    close(f); % close figure after saving

end
