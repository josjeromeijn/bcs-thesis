%%%%%%% Establish good/bad environment growthrate per species. %%%%%%%%
step = 51;

% Load table with E_Y_Good & E_Y_Bad per species
Table_E_Y_Good_Bad = readtable('LOCATION\good_bad.txt');

clear Good_environment Bad_environment
% Create output matrix Good_environment and Bad_environment
Good_environment = array2table(zeros(size(Table_E_Y_Good_Bad,1),8),'VariableNames',{'Species','lambda','R0','GT','Lm','Lp','Lb','Rb'});
Good_environment.Species = Table_E_Y_Good_Bad.id;

Bad_environment = array2table(zeros(size(Table_E_Y_Good_Bad,1),8),'VariableNames',{'Species','lambda','R0','GT','Lm','Lp','Lb','Rb'});
Bad_environment.Species = Table_E_Y_Good_Bad.id;

% Loop over all species
for i = [1:size(Table_E_Y_Good_Bad,1)]
    id = Table_E_Y_Good_Bad.id(i);
    E_Ygood = Table_E_Y_Good_Bad.E_Y_Good(i);
    E_Ybad = Table_E_Y_Good_Bad.E_Y_Bad(i);
    
    %Divide between E_Y = [0.3,0.8] & E_Y = [0.5,1]
    clear All_Species
    if Table_E_Y_Good_Bad.E_Ymin(i) == 0.5        
        All_Species = All_Species_0_5;
    else
        All_Species = All_Species_0_3;
    end
    Lm = lm_table.Lm(i);
    
    %Find file of id
    file_name = 'LOCATION\\id' + string(id) + '.txt';
    single_id_table = readtable(file_name);
    
    %Find lowval (Lb) and Lp of id
    Lb = single_id_table.size_t(1);
    Lp = single_id_table.size_t(end);
    
    %Find rb
    Rb = vB_table.Var2(i);
    
    %Get lambda, R0 and GT of good & bad environment
    [lambda,R0,GT] = get_lambda_R0_GT(All_Species,E_Ygood,i);
    %Put in table
    Good_environment(i,:) = num2cell([id,lambda,R0,GT,Lm,Lp,Lb,Rb]);

    [lambda,R0,GT] = get_lambda_R0_GT(All_Species,E_Ybad,i);
    %Put in table
    Bad_environment(i,:) = num2cell([id,lambda,R0,GT,Lm,Lp,Lb,Rb]);

end

% Remove populations from table (keep 3.1 & 7.3)
indices = [4,5,9,10];
Good_environment(indices,:)=[];
Bad_environment(indices,:)=[];

% Change id to species 
names = readtable('LOCATION\ott_names.txt');
Good_environment.Species = names.OTT_name;
Bad_environment.Species = names.OTT_name;

% Save Good_environment & Bad_environment
writetable(sortrows(Good_environment), 'LOCATION\Good_environment.xlsx')
writetable(sortrows(Bad_environment), 'LOCATION\Bad_environment.xlsx')

