%% Data import

clearvars
close all
clc

df1 = readtable('data_part_1.csv', 'ReadVariableNames',true);
df2 = readtable('data_part_2.csv', 'ReadVariableNames',true);

%%

df1_names = df1.Properties.VariableNames;
df2_names = df2.Properties.VariableNames;


df2_not_in_df1 = setdiff(df2_names, df1_names);
df1_not_in_df2 = setdiff(df1_names, df2_names);

for i = 1:length(df2_not_in_df1)
    missing = df2_not_in_df1{i};
    df1.(missing) = NaN(height(df1), 1);
end


df1_names = df1.Properties.VariableNames;
df2_names = df2.Properties.VariableNames;

df1 = df1(:, sort(df1_names));
df2 = df2(:, sort(df2_names));

df_complete = [df1;df2];