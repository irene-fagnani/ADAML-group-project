%% Data import and combination
clearvars
close all
clc

df1_orig = readtable('data_part_1.csv', 'ReadVariableNames',true);
df2_orig = readtable('data_part_2.csv', 'ReadVariableNames',true);
df1 = df1_orig;
df2 = df2_orig;

% Checking if the column names match
df1_names = df1.Properties.VariableNames;
df2_names = df2.Properties.VariableNames;

df2_not_in_df1 = setdiff(df2_names, df1_names);
df1_not_in_df2 = setdiff(df1_names, df2_names);

for i = 1:length(df2_not_in_df1)
    missing = df2_not_in_df1{i};
    df1.(missing) = NaN(height(df1), 1);
end

df_names = df2.Properties.VariableNames; % Because df2 has all names and they're already in correct order.

df1_complete = df1(:, df_names);
df2_complete = df2(:, df_names);

df_complete = [df1_complete;df2_complete];
df_complete.Var1 = [];

M_complete = table2array(unique(df_complete, 'rows'));
M_traits = M_complete(:,1:37);
M_wl = M_complete(:,38:end);

%% Plotting of trait observations in all data rows.
traits = M_complete(:,1:37);

observations = ones(height(traits), width(traits));

% Loop through each cell in the table
for row = 1:height(traits)
    for col = 1:width(traits)
        cell_value = traits(row, col);
        % Check if the cell is empty (either NaN or empty string)
        if (isnumeric(cell_value) && isnan(cell_value)) || (ischar(cell_value) && isempty(cell_value))
            observations(row, col) = 0;
        end
    end
end

y = sum(observations);
x = 1:1:37;

figure;
bar(x,y)
title('Number of observations per trait variable');
xlabel('Trait variable index');
ylabel('Number of Obsevations');

for i = 1:length(y)
    text(x(i), y(i), num2str(y(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

% Adjust x-axis limits to make sure all text is visible
xlim([0 length(y) + 1]);

%% Plotting all wavelengths for visualization purposes
% Missing points included as NaN values

df_wl = M_complete(:, 38:end);
df_wl_1 = df_wl(:,1:951);
df_wl_2 = df_wl(:,952:1321);
df_wl_3 = df_wl(:,1322:1721);
df_wl = [NaN(length(df_wl), 399),df_wl_1,...
           NaN(length(df_wl),80), df_wl_2,...
           NaN(length(df_wl),250), df_wl_3,...
           NaN(length(df_wl),51)];

x_vect = 1:width(df_wl);
mean_vect = mean(df_wl);
min_vect = min(df_wl);
max_vect = max(df_wl);

figure;
hold on;
plot(mean_vect,'DisplayName','mean');
plot(min_vect, 'DisplayName','minimum');
plot(max_vect, 'DisplayName','maximum')
hold off;

title('Wavelength data')
xlim([0, 2501])
xlabel('Wavelength');
ylabel('Reflectance');
legend('show');
grid on;

%% Additional visualizations:

% Boxplot of trait vars
figure;
boxplot(M_traits)

% Corrplot is not possible for all response variables; there are too many
%% Matrix creations
commonColumns = intersect(df1_orig.Properties.VariableNames, df2_orig.Properties.VariableNames);

df1_common = df1_orig(:, commonColumns);
df2_common = df2_orig(:, commonColumns);

df_tot = [df1_common; df2_common];
df_totM = table2array([df1_common; df2_common]);

Matrices = cell(20, 4);
name_list=["Anth","Boron","C","Ca","Car","Cellulose","Chl","Copper","EWT","Fiber","LAI","Lignin","LMA","Magnesium","Manganese","N","NSC","Phosphorus","Potassium","Sulfur"];
for i = 1:20
    nans = isnan(df_totM(:,i));
    idx=find(nans==0);
    idx_miss=find(nans==1);
    
    X = df_totM(nans==0, 22:end);
    X_nans = df_totM(nans, 22:end);
    Y = df_totM(nans==0, i);
    
    % Scaling and centering
    [X_bar,mu,sigma] = zscore(X);
    X_nans = normalize(X_nans,'center',mu,'scale',sigma);
    
    Matrices{i, 1} = name_list(i);
    Matrices{i, 2} = X; % X with found Y
    Matrices{i, 3} = Y; % Found Y
    Matrices{i, 4} = X_nans; % Scaled and Centered missing y X-values
    Matrices{i, 6} = idx; % Indexes of found values
    Matrices{i, 7} = idx_miss; %indexes of missing values
end

clearvars -except Matrices name_list
close all
clc

%% Visualizing of preprocessed data

figure;
sgtitle("Box-plot of pre-processed data.")
for ii=1:size(Matrices,1)
    subplot(2,10,ii);
    boxplot(Matrices{ii,3});
    xlabel(Matrices{ii,1});
end

%% Bar chart of trait observations
n = length(Matrices);
y = ones(n,1);

% Loop through each cell in the table
for trait = 1:n
    y(trait) = length(Matrices{trait,2});
end

x = 1:20;

figure;
bar(x,y)
title('Number of observations per trait variable', 'FontSize', 30);
xlabel('Trait name', 'FontSize', 25);
ylabel('Number of Obsevations', 'FontSize', 25);
set(gca, 'FontSize', 20);


xticks(1:length(name_list));
xticklabels(name_list);
xtickangle(45);  

for i = 1:length(y)
    text(x(i), y(i), num2str(y(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

% Adjust x-axis limits to make sure all text is visible
xlim([0 length(y) + 1]);

%% PLS/PCR Model and k-fold Cross-Validation
clc
nLV = 50;
k = 5;
bestModels = cell(20, 4); % Store best model results for each trait.

for kk=2
    % Into predictor variables and predicted variable
    X0 = Matrices{kk,2};
    Y0 = Matrices{kk,3};

    % Making division between calibration and validation
    nobs    = length(Y0);
    
    % We save 30% of the data for validation
    part    = cvpartition(nobs, 'HoldOut', 0.3);
    idxCal  = training(part);
    idxTest  = test(part);
    
    X = X0(idxCal, :);
    Y = Y0(idxCal);

    cv = cvpartition(size(Y, 1), 'KFold', k);
    disp(['Training model for trait ', num2str(kk)])
    nCols = size(X0,2)+1;

    bPCR = zeros(nCols, 1);
    PRESSPCR = zeros(k, nLV);
    RMSEPCR = zeros(k, nLV);
    Q2PCR = zeros(k, nLV);
         
    PRESSPLS = zeros(k, nLV);
    RMSEPLS = zeros(k, nLV);
    Q2PLS = zeros(k, nLV);

    for i = 1:k
        trainIdx = training(cv, i);
        testIdx  = test(cv, i);

        X_cal  = X(trainIdx, :);
        Y_cal = Y(trainIdx, :);

        X_val  = X(testIdx, :);
        Y_val  = Y(testIdx, :);

        % Center and Scale
        [XCal, mu, sigma] = zscore(X_cal);
        XVal = normalize(X_val, "center", mu, "scale", sigma);
        YCal = Y_cal - mean(Y_cal);
        YVal = Y_val - mean(Y_cal);

        [P, T, latent] = pca(XCal, 'Centered', false, 'Economy', false);

        TSS = sum((YCal - mean(YCal)).^2);

        [rows, ~] = size(XVal);
        
        
        for j = 1:nLV
            disp(['Trait ', num2str(kk), ', fold ', num2str(i), ',  part ', num2str(j)])
            bPCR(1:end-1)        = P(:,1:j) * regress(YCal, T(:,1:j));
            bPCR(end)            = mean(YCal) - mean(XCal) * bPCR(1:end-1);
            YPredPCR             = [ones(rows,1) XVal] * bPCR;
            PRESSPCR(i,j)        = sum((YPredPCR - YVal).^2);
            RMSEPCR(i,j)         = rmse(YPredPCR, YVal);
            Q2PCR(i,j)           = 1 - PRESSPCR(i,j)/TSS;

            [~, ~, ~, ~, bPLS]   = plsregress(XCal, YCal, j);
            YPredPLS             = [ones(rows,1) XVal] * bPLS;
            PRESSPLS(i,j)        = sum((YPredPLS - YVal).^2);
            RMSEPLS(i,j)         = rmse(YPredPLS, YVal);
            Q2PLS(i,j)           = 1 - PRESSPLS(i,j)/TSS;
        end
    end
    disp(['Training for trait ', num2str(kk), ' completed.'])
    Q2_CV_PLS = mean(Q2PLS);
    Q2_CV_PCR = mean(Q2PCR);
    
    RMSE_CV_PLS = mean(RMSEPLS);
    RMSE_CV_PCR = mean(RMSEPCR);
    
    % LV number selection
    Q2_max_PLS = max(Q2_CV_PLS);
    Q2_mod_PLS = Q2_max_PLS * 0.05;
    Q2_upper_lim_PLS = Q2_max_PLS - Q2_mod_PLS;

    RMSE_min_PLS = min(RMSE_CV_PLS);
    RMSE_mod_PLS = RMSE_min_PLS * 0.05;
    RMSE_lower_lim_PLS = RMSE_min_PLS + RMSE_mod_PLS;

    Opt_noLV = min(min(find(RMSE_CV_PLS < RMSE_lower_lim_PLS)), min(find(Q2_CV_PLS > Q2_upper_lim_PLS)));

    % Component number selection
    Q2_max_PCR = max(Q2_CV_PCR);
    Q2_mod_PCR = Q2_max_PCR * 0.05;
    Q2_upper_lim_PCR = Q2_max_PCR - Q2_mod_PCR;

    RMSE_min_PCR = min(RMSE_CV_PCR);
    RMSE_mod_PCR = RMSE_min_PCR * 0.05;
    RMSE_lower_lim_PCR = RMSE_min_PCR + RMSE_mod_PCR;

    Opt_noComp = min(min(find(RMSE_CV_PCR < RMSE_lower_lim_PCR)), min(find(Q2_CV_PCR > Q2_upper_lim_PCR)));

    if Q2_CV_PLS(Opt_noLV) > Q2_CV_PCR(Opt_noComp)

        if RMSE_CV_PLS(Opt_noLV) < RMSE_CV_PCR(Opt_noComp)
            bestModel = 'PLS';
            bestRMSE = RMSE_CV_PLS(Opt_noLV);
            bestQ2 = Q2_CV_PLS(Opt_noLV);

        elseif abs(RMSE_CV_PLS(Opt_noLV) - RMSE_CV_PCR(Opt_noComp)) > 0.05
            bestModel = 'PCR';
            bestRMSE = RMSE_CV_PCR(Opt_noComp);
            bestQ2 = Q2_CV_PCR(Opt_noComp);

        end

    else
        if RMSE_CV_PLS(Opt_noLV) > RMSE_CV_PCR(Opt_noComp)
            bestModel = 'PCR';
            bestRMSE = RMSE_CV_PCR(Opt_noComp);
            bestQ2 = Q2_CV_PCR(Opt_noComp);

        elseif abs(RMSE_CV_PLS(Opt_noLV) - RMSE_CV_PCR(Opt_noComp)) > 0.05
            bestModel = 'PLS';
            bestRMSE = RMSE_CV_PLS(Opt_noLV);
            bestQ2 = Q2_CV_PLS(Opt_noLV);

        end
    end
    
    % Final predictions
    XPred = Matrices{kk, 4}; % X matrix with missing trait values
    X0 = zscore(X0); % Scale X0 to match the previously scaled XPred

    if strcmp(bestModel,'PLS')
        % Predict using PLS
        [~, ~, ~, ~, bPLS] = plsregress(X0, Y0, Opt_noLV);
        YPred = [ones(size(XPred, 1), 1), XPred] * bPLS; % Predict using PLS.
    else
         % Predict using PCR.
        [P, T, ~] = pca(X0, 'Centered', false, 'Economy', false);
        bPCR = P(:, 1:Opt_noComp) * regress(Y0 - mean(Y0), T(:, 1:Opt_noComp));
        bPCR = [mean(Y0) - mean(X0) * bPCR; bPCR];
        YPred = [ones(size(XPred, 1), 1), XPred] * bPCR;
    end
    
    Matrices{kk, 5} = YPred;
    
    % Store the selected model information
    bestModels{kk, 1} = name_list(kk);
    bestModels{kk, 2} = bestModel;
    bestModels{kk, 3} = bestRMSE;
    bestModels{kk, 4} = bestQ2;

    figure;
    hold on;
    plot(Q2_CV_PLS);
    plot(Q2_CV_PCR);
    plot(Opt_noLV, Q2_CV_PLS(Opt_noLV), 'b.', 'MarkerSize', 40)
    plot(Opt_noComp, Q2_CV_PCR(Opt_noComp), 'r.', 'MarkerSize', 40)
    title('Q2 comparison', 'FontSize', 30)
    xlabel("No of LVs/Components in the model.", 'FontSize', 25)
    ylabel("Q^2_{CV}", 'FontSize', 25)
    legend(["PLS"; "PCR";""], 'FontSize', 20)
    set(gca, 'FontSize', 20);

    figure;
    hold on;
    plot(RMSE_CV_PLS);
    plot(RMSE_CV_PCR);
    plot(Opt_noLV, RMSE_CV_PLS(Opt_noLV), 'b.', 'MarkerSize', 40)
    plot(Opt_noComp, RMSE_CV_PCR(Opt_noComp), 'r.', 'MarkerSize', 40)
    title('RMSE comparison', 'FontSize', 30)
    xlabel("No of LVs/Components in the model.", 'FontSize', 25)
    ylabel("RMSE_{CV}", 'FontSize', 25)
    legend(["PLS"; "PCR"], 'FontSize', 20);
    set(gca, 'FontSize', 20);
end

% Display the selected best models for each trait
disp('Best models for each trait:');
disp(cell2table(bestModels, 'VariableNames', {'Trait', 'BestModel', 'BestRMSE', 'BestQ2'}));

%% Create the complete trait matrix
trait_table = zeros(13295, 20);
for i = 1:length(Matrices)
    found_idx = Matrices{i, 6};
    missing_idx = Matrices{i, 7};
    trait_table(found_idx,i) = Matrices{i, 3};
    trait_table(missing_idx,i) = Matrices{i, 5};
end
filename = "trait_matrix.mat";
save(filename, 'trait_table')

%% PCA on the complete trait matrix
clearvars
close all -except Matrices
clc


name_list=["Anth","Boron","C","Ca","Car","Cellulose","Chl","Copper","EWT","Fiber","LAI","Lignin","LMA","Magnesium","Manganese","N","NSC","Phosphorus","Potassium","Sulfur"];

trait_matrix_data = load('trait_matrix.mat');
trait_matrix = trait_matrix_data.trait_table;

figure;
boxplot(trait_matrix)
title('Boxplot of trait data', 'FontSize', 30)
xticks(1:length(name_list));
xticklabels(name_list);
xtickangle(45);
set(gca, 'FontSize', 20);

% Rescale and visualize again
trait_matrix = zscore(trait_matrix);

figure;
boxplot(trait_matrix)
title('Boxplot of scaled trait data', 'FontSize', 30)
xticks(1:length(name_list));
xticklabels(name_list);
xtickangle(45);
set(gca, 'FontSize', 20);

% PCA using MatLab's pca()-function
[Loadings, Scores, EigenVals, T2, Explained, mu] = pca(trait_matrix);

% Plotting cumulative explained variance
figure;
pareto(Explained)
title('Cumulative explained variance', 'FontSize', 30)
xlabel('Principal component', 'FontSize', 25);
ylabel('% of variance explained', 'FontSize', 25)
set(gca, 'FontSize', 20);

% Appears that 9 principal components are required to explain 95% of cumulative variance
name_list=["Anth","Boron","C","Ca","Car","Cellulose","Chl","Copper","EWT","Fiber","LAI","Lignin","LMA","Magnesium","Manganese","N","NSC","Phosphorus","Potassium","Sulfur"];

% biplot
figure;
biplot(Loadings(:, 1:2), 'scores', Scores(:, 1:2), 'varlabels', name_list)
title('Biplot of traits', 'FontSize', 30)
set(gca, 'FontSize', 20);

% Visualizing Loading values

% Component 1 Loadings
figure;
bar(Loadings(:,1));
title('Component 1 Loadings', 'FontSize', 30);
xlabel('Trait', 'FontSize', 25);
ylabel('Loading value', 'FontSize', 25);
set(gca, 'FontSize', 20);

xticks(1:length(name_list));
xticklabels(name_list);
xtickangle(45);                           

% Component 2 Loadings
figure;
bar(Loadings(:,2));
title('Component 2 Loadings', 'FontSize', 30);
xlabel('Trait', 'FontSize', 25);
ylabel('Loading value', 'FontSize', 25);
set(gca, 'FontSize', 20);

xticks(1:length(name_list));              
xticklabels(name_list);                   
xtickangle(45);      
%% Get the coefficients beta
%noPCsPCR =  Opt_noComp;
%noPCsPLS =   Opt_noLV;

noPCsPCR =  10;
noPCsPLS =   15;

% considering XCal all the Matrices{kk,2} and XVal the Matrices{kk,4}

for kk=1:20
XCal=Matrices{kk,2};
YCal=Matrices{kk,3};
[P, T, latent] = pca(XCal, 'Centered', false, 'Economy', false);
% Re-calibrate the models with the combined cal-val partition
bPCR = [];
bPCR = P(:,1:noPCsPCR) * regress(YCal, T(:,1:noPCsPCR)); 
bPCR = [mean(YCal) - mean(XCal) * bPCR; bPCR]; % Add intercept.
     
[~, ~, ~, ~, bPLS] = plsregress(XCal, YCal, noPCsPLS);

figure;
betas = [bPLS(2:end), bPCR(2:end)];
bar(betas);
legend(["PLS Regression Coefficients", "PCR Regression Coefficients"]);
end

%% Residuals
% resid
residPLS = abs(YVal - YTestPredPLS);
residPCR = abs(YVal - YTestPredPCR);

% subplot(1,2,1);
% scatter(model(1).YTest, residPLS, 'filled');
% xlabel("Age [normalized]");
% ylabel("Estimate Age error");
% title("PLS Estimation Residuals");

% subplot(1,2,2);
% scatter(model(1).YTest, residPCR, 'filled');
% xlabel("Age [normalized]");
% ylabel("Estimate Age error");
% title("PCR Estimation Residuals");
