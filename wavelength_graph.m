%% graphs (wavelength part)
clear("df1")
clear("df2")

df1 = readtable('data_part_1.csv', 'ReadVariableNames',true);
df2 = readtable('data_part_2.csv', 'ReadVariableNames',true);

df1_wl=df1(:,22:1742);
df2_wl=df2(:,39:1759);
%%
df_wl=[df1_wl;df2_wl];
df_wl=table2array(df_wl);
%% add NaN columns where water absorption happens (1351–1430, 1801–2050 and 2451–2501 nm)

%% first bandwith
col_nan=NaN(13295,1);

for ii=1:80 % the 952 column should be the first NaN column
    df_wl = [df_wl(:, 1:950+ii) col_nan df_wl(:, 951+ii:end)];
end

size(df_wl,2)

%% second bandwith
pos=1322+79-1;

for ii=1:250 
    df_wl = [df_wl(:, 1:pos+ii) col_nan df_wl(:, pos+ii+1:end)];
end


%%

% substitute min and max with 99% quantiles
mean_vect=(mean(df_wl));
min_vect=prctile(df_wl,1);
max_vect=prctile(df_wl,99);

%%
x_vect=linspace(400,2450,2051);
hold on;
plot(x_vect,mean_vect,'DisplayName','mean');
plot(x_vect, min_vect, 'DisplayName','1% percentile');
plot(x_vect, max_vect, 'DisplayName','99% percentile')
hold off;

xlabel('wavelength');
ylabel('Reflectance');
legend('show');
grid on;

%% boxplot