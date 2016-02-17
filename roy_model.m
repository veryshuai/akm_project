%Simulate a Roy model in wages

rng(80083)

N_workers=2000;
N_firms=200;
N_w_types=10;
N_f_types=10;
N_years=1000;
match_effect_scale=-0.2; %strength of the match effect
util_disp=2; %dispersion of worker-firm specific utility draw
unemp_prob=0.1; %probability of being knocked to unemp

%Get permanent firm and worker types
firm_type_list = [randi(N_f_types,N_firms,1); -inf]; %last spot is for unemployed
worker_type_list = randi(N_w_types,N_workers,1);
match_util_effect = [normrnd(0,util_disp,N_workers,N_firms), ones(N_workers,1) * -inf]; %extra column is for unemployed

%Worker placement vector
workers_firm = ones(N_workers,1) * (N_firms + 1); %unemployed
l_wage = zeros(N_workers,1); %initial wage vector
util = ones(N_workers,1) * -inf; %initial utility

%simulate 1000 years
workers_firm_by_year = zeros(N_workers,N_years);
workers_wage_by_year = zeros(N_workers,N_years);
for year=1:N_years
    
    %Unemployment for some
    unemployed = rand(N_workers,1) < unemp_prob;
    workers_firm(unemployed) = N_firms + 1;
    util(unemployed) = -inf;
    l_wage(unemployed) = -inf;
    
    %Each worker gets an offer from a random firm
    offer_firm = randi(N_firms,N_workers,1);
    
    %Wage from new firm
    l_wage_offer = worker_type_list + firm_type_list(offer_firm) + match_effect_scale *  worker_type_list .* firm_type_list(offer_firm);
    util_offer = match_util_effect(sub2ind(size(match_util_effect),[1:length(offer_firm)]',offer_firm)) + l_wage_offer; %choose the relevant match specific utility and wage
    
    %Accept?  If so update all the states
    accepters = util_offer > util;
    workers_firm(accepters) = offer_firm(accepters);
    util(accepters) = util_offer(accepters);
    l_wage(accepters) = l_wage_offer(accepters);
    
    %Record state
    workers_firm_by_year(:,year) = workers_firm;
    workers_wage_by_year(:,year) = l_wage;
    
end

%Perform Card Test

%Calculate mean coworker wages by year for last 8 years
tot_wages_by_firm = zeros(N_firms,8);
emp_by_firm = zeros(N_firms,8);
cow_wages_by_worker = zeros(N_workers,8);
for year=1:8
    result_year_index = N_years - year + 1; %get the index for year from results
    tot_wages_by_firm(:,year) = accumarray([workers_firm_by_year(:,result_year_index) ones(N_workers,1)],workers_wage_by_year(:,result_year_index));
    emp_by_firm(:,year) = accumarray([workers_firm_by_year(:,N_years-year+1) ones(N_workers,1)],ones(N_workers,1));
    cow_wages_by_worker(:,year) = (tot_wages_by_firm(workers_firm_by_year(:,result_year_index),year) - workers_wage_by_year(:,result_year_index)) ./ (emp_by_firm(workers_firm_by_year(:,result_year_index),year) - 1); 
end

%Put firms into four groups based on coworker wages
firm_groups = ones(N_workers,8);
for year=1:8
    cutoffs = quantile(cow_wages_by_worker(:,year),3);
    for k=1:3
        above_quant = cow_wages_by_worker(:,year) > cutoffs(k); 
        firm_groups(above_quant,year) = firm_groups(above_quant,year) + 1;
    end
end

%Find moves between 1 and 4 and 4 and 1
movin_up = firm_groups(:,1:end-1) == 1 & firm_groups(:,2:end) == 4;
movin_up_source_wage = movin_up .* workers_wage_by_year(:,end-8:end-2);
movin_up_dest_wage = movin_up .* workers_wage_by_year(:,end-7:end-1);
movin_up_change_in_average = sum(sum(movin_up_dest_wage)) / sum(sum(movin_up)) - sum(sum(movin_up_source_wage)) / sum(sum(movin_up));
movin_down = firm_groups(:,1:end-1) == 4 & firm_groups(:,2:end) == 1;
movin_down_source_wage = movin_down .* workers_wage_by_year(:,end-8:end-2);
movin_down_dest_wage = movin_down .* workers_wage_by_year(:,end-7:end-1);
movin_down_change_in_average = sum(sum(movin_down_dest_wage)) / sum(sum(movin_down)) - sum(sum(movin_down_source_wage)) / sum(sum(movin_down));

%Display results
display(movin_up_change_in_average);
display(movin_down_change_in_average);

%Do AKM
wid = repmat((1:size(workers_firm))',1,8); %create a worker id variables
l_wage_reg = workers_wage_by_year(:,end-7:end); %only required years
fid = workers_firm_by_year(:,end-7:end); %only the years we want
data = mat2dataset([l_wage_reg(:),wid(:),fid(:)],'VarNames',{'L_wage','WID','FID'}); %put data in a matlab statistics readable dataframe
data.WID = nominal(data.WID); %make this a category
data.FID = nominal(data.FID); %make this a category
model = fitlme(data,'L_wage ~ WID + FID');
[beta,betanames,stats] = fixedEffects(model);

%How did it do?
AKM_est_workers = beta(2:N_workers);
actual_workers = worker_type_list(2:end);
display(corr(actual_workers,AKM_est_workers))
AKM_est_firms = beta(N_workers+1:end);
actual_firms = firm_type_list(2:end-1);
display(corr(AKM_est_firms,actual_firms));


