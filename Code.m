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
    [X_bar,mu,sigma]=zscore(X);
    X_new=normalize(X,'center',mu,'scale',sigma);
    Y = df_totM(nans==0, i);
    X_nans = df_totM(nans, 22:end);

    Matrices{i, 1} = name_list(i);
    Matrices{i, 2} = zscore(X_new);
    Matrices{i, 3} = zscore(Y);
    Matrices{i, 4} = zscore(X_nans);
    Matrices{i, 5}=idx; 
    Matrices{i, 6}=idx_miss; %indexes of missing values
end


%% Visualizing of preprocessed data

figure;
sgtitle("Box-plot of pre-processed data.")
for ii=1:size(Matrices,1)
    subplot(2,10,ii);
    boxplot(Matrices{ii,3});
    xlabel(Matrices{ii,1});
    ylim([-2 14]);
end



%%Ekhtear

%% PLS/PCA Model and Cross-Validation

for kk=1:size(Matrices,1)
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

k = 5;
cv = cvpartition(size(Y, 1), 'KFold', k);


for i = 1:cv.NumTestSets

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

    for j = 1:19

        bPCR = [];
        bPCR            = P(:,1:j) * regress(YCal, T(:,1:j));
        bPCR            = [mean(YCal) - mean(XCal) * bPCR; bPCR];
        YPredPCR        = [ones(rows,1) XVal] * bPCR;
        PRESSPCR(i,j)   = sum((YPredPCR - YVal).^2);
        RMSEPCR(i,j)    = rmse(YPredPCR, YVal);
        Q2PCR(i,j)      = 1 - PRESSPCR(i,j)/TSS;

        [~, ~, ~, ~, bPLS] = plsregress(XCal, YCal, j);
        YPredPLS        = [ones(rows,1) XVal] * bPLS;
        PRESSPLS(i,j)   = sum((YPredPLS - YVal).^2);
        Q2PLS(i,j)      = 1 - PRESSPLS(i,j)/TSS;
        RMSEPLS(i,j)    = rmse(YPredPLS, YVal);

    end

end
end

%%
Q2_CV_PCR = mean(Q2PCR);
Q2_CV_PLS = mean(Q2PLS);

PRESS_CV_PLS = mean(PRESSPLS);
PRESS_CV_PCR = mean(PRESSPCR);

RMSE_CV_PLS = mean(RMSEPLS);
RMSE_CV_PCR = mean(RMSEPCR);

figure;
plot(Q2_CV_PLS);
hold on;
plot(Q2_CV_PCR);
xlabel("No LVs in the model.")
ylabel("Q^2_{CV}")
legend(["PLS"; "PCR"])

figure;
plot(RMSE_CV_PLS);
hold on;
plot(RMSE_CV_PCR);
xlabel("No LVs in the model.")
ylabel("RMSE_{CV}")
legend(["PLS"; "PCR"]);
















%% Predict Missing Trait Values for New Data
% normalize to scale new data using the saved mean and std
% X_new = normalize(M_wl, 'center', X_train_mean, 'scale', X_train_std);
