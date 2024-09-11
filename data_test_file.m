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
df_complete.Var1 = [];
df_complete = table2array(unique(df_complete, 'rows'));

%%
traits = df_complete(:,1:37);

empty_cells = ones(height(traits), width(traits));

% Loop through each cell in the table
for row = 1:height(traits)
    for col = 1:width(traits)
        cell_value = traits(row, col);
        % Check if the cell is empty (either NaN or empty string)
        if (isnumeric(cell_value) && isnan(cell_value)) || (ischar(cell_value) && isempty(cell_value))
            empty_cells(row, col) = 0;
        end
    end
end

y = sum(empty_cells)
x = 1:1:37

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

%% graphs
clear("df1")
clear("df2")

df1 = readtable('data_part_1.csv', 'ReadVariableNames',true);
df2 = readtable('data_part_2.csv', 'ReadVariableNames',true);

df1_wl=df1(:,22:1742);
df2_wl=df2(:,39:1759);

df_wl=[df1_wl;df2_wl];
df_wl=table2array(df_wl);
%% add NaN columns where water absorption happens (1351-1430; 1801-2023;2451-2501)
col_nan=NaN(13295,1);

for ii=1:80 % the 952 column should be the first NaN column
    df_wl = [df_wl(:, 1:950+ii) col_nan df_wl(:, 951+ii+1:end)];
end

pos=80+(1801-401);

for ii=1:223 
    df_wl = [df_wl(:, 1:pos+ii) col_nan df_wl(:, pos+ii+1:end)];
end


%%
x_vect=linspace(400,2450, 1721);

mean_vect=(mean(df_wl));
min_vect=(min(df_wl));
max_vect=(max(df_wl));

%%
hold on;
plot(mean_vect,'DisplayName','mean');
plot(min_vect, 'DisplayName','minimum');
plot(max_vect, 'DisplayName','maximum')
hold off;

xlabel('wavelength');
ylabel('Reflectance');
legend('show');
grid on;

%% boxplot




