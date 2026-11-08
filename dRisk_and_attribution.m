models = {'ACCESS-ESM1-5','BCC-CSM2-MR','CanESM5','CESM2-WACCM','CMCC-CM2-SR5',...
    'CMCC-ESM2','EC-Earth3-Veg','EC-Earth3-Veg-LR','INM-CM4-8','INM-CM5-0',...
    'IPSL-CM6A-LR','KIOST-ESM','MPI-ESM1-2-HR','MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','TaiESM1'};
scenarios = {'SSP245', 'SSP585'};
vars = {'TAS', 'PR'};
%
years = 1980:2100;
n_years = numel(years);

month_edges = datetime(1980, 1, 1) + calmonths(0:12*n_years);
days_per_month = days(diff(month_edges));

seconds_per_month = reshape(days_per_month * 86400,1, 1, 12, n_years);
data_struct = nan(180, 360, n_years, ...
                  length(scenarios), length(vars), length(models));

for s = 1:length(scenarios)
    for v = 1:length(vars)
        for m = 1:length(models)

            file_name = sprintf( ...
                'G:/data and documents/CMIP6TP/%s%smonth_%s_1980-2100.mat', ...
                vars{v}, scenarios{s}, models{m});

            loaded_data = load(file_name);
            fields = fieldnames(loaded_data);
            temp_data = loaded_data.(fields{1});
            if v == 1
                annual_data = mean(temp_data, 3, 'omitnan') - 273.15;
            else
                monthly_precip = temp_data .* seconds_per_month;
                annual_data = sum(monthly_precip, 3);
            end

            data_struct(:, :, :, s, v, m) = ...
                reshape(annual_data, 180, 360, n_years);
        end
    end
end
%%
load crowther_stem.mat %tree number distribution
future_mat=squeeze(data_struct(:,:,43:121,1,1,:));%2023-2100 climate
future_map=squeeze(data_struct(:,:,43:121,1,2,:));

data = readtable('melted_data_hot.csv');
tree_mask=crowther_resampled;tree_mask(isnan(tree_mask))=0;
value_h245_ci = zeros(180, 360,79,17);
        for i = 1:180
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:,:)); % (°C)
                    prec = squeeze(future_map(i, j,:,:)); % (mm)
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1); % -20 ~ 30°C
                    col_idx = max(min(ceil(prec /10),500),1); % 0 ~ 5000mm
                    for iyr=1:79
                    for imod=1:17
                    idx = find(data.Row == row_idx(iyr,imod) & data.Column == col_idx(iyr,imod));
                    if ~isempty(idx)
                        value_h245_ci(i, j,iyr,imod) = data.Value(idx);
                    end
                    end
                    end
                    end
                end
            end
        end


% SSP585

future_mat=squeeze(data_struct(:,:,43:121,2,1,:));%2023-2100 climate
future_map=squeeze(data_struct(:,:,43:121,2,2,:));

value_h585_ci = zeros(180, 360,79,17);
        for i = 1:180
            %disp(i)
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:,:)); 
                    prec = squeeze(future_map(i, j,:,:)); 
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1); 
                    col_idx = max(min(ceil(prec /10),500),1);
                    for iyr=1:79
                    for imod=1:17
                    idx = find(data.Row == row_idx(iyr,imod) & data.Column == col_idx(iyr,imod));
                    if ~isempty(idx)
                        value_h585_ci(i, j,iyr,imod) = data.Value(idx);
                    end
                    end
                    end
                    end
                end
            end
        end

save SSP_hot_17model.mat value_h245_ci value_h585_ci
clear value_h*

%
future_mat=squeeze(data_struct(:,:,43:121,1,1,:));%2023-2100 climate
future_map=squeeze(data_struct(:,:,43:121,1,2,:));

data = readtable('melted_data_cold.csv');
tree_mask=crowther_resampled;tree_mask(isnan(tree_mask))=0;
value_c245_ci = zeros(180, 360,79,17);
        for i = 1:180
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:,:)); 
                    prec = squeeze(future_map(i, j,:,:)); 
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1); 
                    col_idx = max(min(ceil(prec /10),500),1); 
                    for iyr=1:79
                    for imod=1:17
                    idx = find(data.Row == row_idx(iyr,imod) & data.Column == col_idx(iyr,imod));
                    if ~isempty(idx)
                        value_c245_ci(i, j,iyr,imod) = data.Value(idx);
                    end
                    end
                    end
                    end
                end
            end
        end


% SSP585

future_mat=squeeze(data_struct(:,:,43:121,2,1,:));%2023-2100 climate
future_map=squeeze(data_struct(:,:,43:121,2,2,:));

value_c585_ci = zeros(180, 360,79,17);
        for i = 1:180
            %disp(i)
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:,:)); 
                    prec = squeeze(future_map(i, j,:,:)); 
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1);
                    col_idx = max(min(ceil(prec /10),500),1); 
                    for iyr=1:79
                    for imod=1:17
                    idx = find(data.Row == row_idx(iyr,imod) & data.Column == col_idx(iyr,imod));
                    if ~isempty(idx)
                        value_c585_ci(i, j,iyr,imod) = data.Value(idx);
                    end
                    end
                    end
                    end
                end
            end
        end

save SSP_cold_17model.mat value_c245_ci value_c585_ci
clear value_c*

%%
value_h245=tsnanmedian(value_h245_ci,4);value_h585=tsnanmedian(value_h585_ci,4);

num_a1=nan(1,79);num_a1_ci=nan(17,79);num_a2=nan(1,79);num_a2_ci=nan(17,79);
for iyr=1:79
tmp=value_h245(:,:,iyr);
num_a1(iyr)=sum(tree_mask(tmp<0),'all');
tmp=value_h585(:,:,iyr);
num_a2(iyr)=sum(tree_mask(tmp<0),'all');
for imod=1:17
tmp=value_h245_ci(:,:,iyr,imod);
num_a1_ci(imod,iyr)=sum(tree_mask(tmp<0),'all');
tmp=value_h585_ci(:,:,iyr,imod);
num_a2_ci(imod,iyr)=sum(tree_mask(tmp<0),'all');
end
end


%%
value_c245=tsnanmedian(value_c245_ci,4);value_c585=tsnanmedian(value_c585_ci,4);

num_a1=nan(1,79);num_a1_ci=nan(17,79);num_a2=nan(1,79);num_a2_ci=nan(17,79);
for iyr=1:79
tmp=value_c245(:,:,iyr);
num_a1(iyr)=sum(tree_mask(tmp<0),'all');
tmp=value_c585(:,:,iyr);
num_a2(iyr)=sum(tree_mask(tmp<0),'all');
for imod=1:17
tmp=value_c245_ci(:,:,iyr,imod);
num_a1_ci(imod,iyr)=sum(tree_mask(tmp<0),'all');
tmp=value_c585_ci(:,:,iyr,imod);
num_a2_ci(imod,iyr)=sum(tree_mask(tmp<0),'all');
end
end
%% ΔRisk
future_mat=squeeze(tsnanmean(data_struct(:,:,92:121,1,1,:),3));%2023-2100 climate
future_map=squeeze(tsnanmean(data_struct(:,:,92:121,1,2,:),3));

data = readtable('melted_data_hot.csv');
tree_mask=crowther_resampled;tree_mask(isnan(tree_mask))=0;
value_h245_ci = nan(180, 360,17);
        for i = 1:180
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:));
                    prec = squeeze(future_map(i, j,:)); 
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1);
                    col_idx = max(min(ceil(prec /10),500),1);
                    for imod=1:17
                    idx = find(data.Row == row_idx(imod) & data.Column == col_idx(imod));
                    if ~isempty(idx)
                        value_h245_ci(i, j,imod) = data.Value(idx);
                    end
                    end
                    end
                end
            end
        end

future_mat=squeeze(tsnanmean(data_struct(:,:,92:121,2,1,:),3));
future_map=squeeze(tsnanmean(data_struct(:,:,92:121,2,2,:),3));

value_h585_ci = nan(180, 360, 17);
        for i = 1:180
            %disp(i)
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:));
                    prec = squeeze(future_map(i, j,:));
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1);
                    col_idx = max(min(ceil(prec /10),500),1);
                    for imod=1:17
                    idx = find(data.Row == row_idx(imod) & data.Column == col_idx(imod));
                    if ~isempty(idx)
                        value_h585_ci(i, j,imod) = data.Value(idx);
                    end
                    end
                    end
                end
            end
        end

save SSP_hot_17model_last30.mat value_h245_ci value_h585_ci


%
future_mat=squeeze(tsnanmean(data_struct(:,:,92:121,1,1,:),3));%2023-2100 climate
future_map=squeeze(tsnanmean(data_struct(:,:,92:121,1,2,:),3));

data = readtable('melted_data_cold.csv');
tree_mask=crowther_resampled;tree_mask(isnan(tree_mask))=0;
value_c245_ci = nan(180, 360,17);
        for i = 1:180
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:)); 
                    prec = squeeze(future_map(i, j,:));
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1);
                    col_idx = max(min(ceil(prec /10),500),1);
                    for imod=1:17
                    idx = find(data.Row == row_idx(imod) & data.Column == col_idx(imod));
                    if ~isempty(idx)
                        value_c245_ci(i, j, imod) = data.Value(idx);
                    end
                    end                    
                    end
                end
            end
        end


% SSP585

future_mat=squeeze(tsnanmean(data_struct(:,:,92:121,2,1,:),3));%2023-2100 climate
future_map=squeeze(tsnanmean(data_struct(:,:,92:121,2,2,:),3));

value_c585_ci = nan(180, 360,17);
        for i = 1:180
            for j = 1:360
                if tree_mask(i, j)
                    temp = squeeze(future_mat(i, j,:));
                    prec = squeeze(future_map(i, j,:)); 
                    if ~isnan(temp(1))&& ~isnan(prec(1))
                    row_idx = max(min(ceil((temp + 20)*10),500),1);
                    col_idx = max(min(ceil(prec /10),500),1); 
                    for imod=1:17
                    idx = find(data.Row == row_idx(imod) & data.Column == col_idx(imod));
                    if ~isempty(idx)
                        value_c585_ci(i, j, imod) = data.Value(idx);
                    end
                    end
                    end
                end
            end
        end

save SSP_cold_17model_last30.mat value_c245_ci value_c585_ci

%% deltavalue

mask=tree_mask;mask(mask==0)=nan;mask(isnan(MAT_new))=nan;
value_gridh(isnan(mask))=nan;
value_gridc(isnan(mask))=nan;
value_gridh245=tsnanmedian(value_h245_ci,3);value_gridh245(isnan(mask))=nan;
value_gridh585=tsnanmedian(value_h585_ci,3);value_gridh585(isnan(mask))=nan;
value_gridc245=tsnanmedian(value_c245_ci,3);value_gridc245(isnan(mask))=nan;
value_gridc585=tsnanmedian(value_c585_ci,3);value_gridc585(isnan(mask))=nan;

dconf_h245=value_gridh245-value_gridh;
dconf_h585=value_gridh585-value_gridh;
dconf_c245=value_gridc245-value_gridc;
dconf_c585=value_gridc585-value_gridc;

dconf_h245(abs(dconf_h245)<0.001)=0;
dconf_h585(abs(dconf_h585)<0.001)=0;
dconf_c245(abs(dconf_c245)<0.001)=0;
dconf_c585(abs(dconf_c585)<0.001)=0;
%% significance-test
tmp=nan(180,360,17);
for i=1:17
    tmp(:,:,i)=value_h585_ci(:,:,i)-value_gridh;
end
tmp(abs(tmp)<0.001)=0;
p_map_h585=nan(180,360);
for i=1:180
    for j=1:360
val=squeeze(tmp(i,j,:));
val=sign(val(~isnan(val)&val~=0));
sign0=sign(dconf_h585(i,j));
p_map_h585(i,j)=numel(find(val==sign0));
    end
end


sig_mask_h585=zeros(180,360);
sig_mask_h585(p_map_h585>12)=1;


%% significance-test;SSP585

tmp=nan(180,360,17);
for i=1:17
    tmp(:,:,i)=value_c585_ci(:,:,i)-value_gridc;
end
tmp(abs(tmp)<0.001)=0;
p_map_c585=nan(180,360);
for i=1:180
    for j=1:360
val=squeeze(tmp(i,j,:));
val=sign(val(~isnan(val)&val~=0));
sign0=sign(dconf_c585(i,j));
p_map_c585(i,j)=numel(find(val==sign0));
    end
end


sig_mask_c585=zeros(180,360);
sig_mask_c585(p_map_c585>12)=1;


%% extended fig5c-f(dominant factors)

data= readtable('melted_data_hot.csv');
value= reshape(table2array(data(:,3)),[500,500]);

dMAT = -1; dMAP = 100;[nMAT, nMAP] = size(value);

S_MAT0 = nan(nMAT, nMAP);S_MAP0 = nan(nMAT, nMAP);

for i = 2:nMAT-1
    S_MAT0(i,:) = (value(i+1,:) - value(i-1,:)) / (2*dMAT);
end

S_MAT0(1,:)     = (value(2,:)   - value(1,:))   / dMAT;
S_MAT0(end,:)   = (value(end,:) - value(end-1,:)) / dMAT;

for j = 2:nMAP-1
    S_MAP0(:,j) = (value(:,j+1) - value(:,j-1)) / (2*dMAP);
end

S_MAP0(:,1)     = (value(:,2)   - value(:,1))   / dMAP;
S_MAP0(:,end)   = (value(:,end) - value(:,end-1)) / dMAP;

tree_mask=crowther_resampled;tree_mask(isnan(tree_mask))=0;
S_MAT=nan(180,360);S_MAP=nan(180,360);%world-map
for i = 1:180
            for j = 1:360
                if tree_mask(i, j)
                    temp = MAT_new(i, j);
                    prec = MAP_new(i, j); 
                    row_idx = max(min(ceil((temp + 20)*10),500),1);
                    col_idx = max(min(ceil(prec /10),500),1);
                    if ~isempty(idx)
                        S_MAT(i,j) = S_MAT0(row_idx,col_idx);
                        S_MAP(i,j) = S_MAP0(row_idx,col_idx);
                    end
                end
            end
end

%dMAT/dMAP contribuion，multi-model median
dconf = dconf_h245;
future_mat=squeeze(tsnanmean(data_struct(:,:,92:121,1,1,:),3));%2023-2100 climate
future_map=squeeze(tsnanmean(data_struct(:,:,92:121,1,2,:),3));
dom_ratio=nan(180,360,17);
for imdl=1:17
dMAT=future_mat(:,:,imdl)-MAT_new;dMAT(dconf_h245==0)=nan;dMAT(isnan(dconf_h245))=nan;
dMAP=future_map(:,:,imdl)-MAP_new;dMAP(dconf_h245==0)=nan;dMAP(isnan(dconf_h245))=nan;
C_MAT = (S_MAT .* dMAT); C_MAP = (S_MAP .* dMAP);

dom_factor = zeros(size(dconf));
dom_factor(C_MAT > (C_MAP+C_MAT)*0.6 ) = 1;
dom_factor(C_MAP > (C_MAP+C_MAT)*0.6 ) = 2;
dom_factor(isnan(C_MAT)|isnan(C_MAP))=nan;
area_MAT = sum(dom_factor(:) == 1, 'omitnan');
area_MAP = sum(dom_factor(:) == 2, 'omitnan');
area_mix = sum(dom_factor(:) == 0, 'omitnan');

total = area_MAT + area_MAP + area_mix;

dom_ratio(:,:,imdl)=C_MAT./(C_MAT+C_MAP);
end
figure();imagesc(tsnanmedian(dom_ratio,3));colorbar;colormap([1,1,1;othercolor('BuOr_8')*0.9]);clim([0 1])
dom_ratio_ci=dom_ratio;
dom_ratio=tsnanmedian(dom_ratio,3);