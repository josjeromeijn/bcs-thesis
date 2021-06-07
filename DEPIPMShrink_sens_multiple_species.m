% Run a sensitivity analysis over multiple organisms 
% Input parameters for sensitivity analysis: MostSensPar  = DEBIPMShrink_MicroOrg_sens(E_Y,MatrixSize,Lb,Lp,Lm,E_Ystdev,rb,kappa,mu)


% These are set parameters for all species
MatrixSize = 200;
E_Ystdev = 0.1;
kappa = 1;
mu = 0.01;

%Read in all Lm of species
lm_table = readtable('LOCATION\\lm.txt');

%Read in all vB growth rates of species
vB_table = readtable('LOCATION\\growthrates.txt');

%Read in all feeding levels of species
E_Y_table = readtable('LOCATION\\good_bad.txt');

%Create list of all species id's
all_id = lm_table.id;

%Create loop over all species id's
for i = [1:length(all_id)]
    id = all_id(i);
    
    %Find E_Y of id
    E_Y = E_Y_table.E_Y_Good(i);
    
    %Find Lm of id
    Lm = lm_table.Lm(i);
    
    %Find file of id
    file_name = 'LOCATION\\id' + string(id) + '.txt';
    single_id_table = readtable(file_name);
    
    %Find Lb and Lp of id
    Lb = single_id_table.size_t(1);
    Lp = single_id_table.size_t(end);
    
    %Find rb of id
    rb = vB_table.Var2(i);
    
    %Do sensitivity analysis
    MostSensPar  = DEBIPMShrink_MicroOrg_sens(E_Y,MatrixSize,Lb,Lp,Lm,E_Ystdev,rb,kappa,mu);
    
    %Save outcomes in table
    Sensitivity_All_Species(i,:) = [id,MostSensPar];
end

%Save Sensitivity_All_Species
writematrix(Sensitivity_All_Species,'LOCATION\\DEB_IPM_sensitivity.txt','Delimiter','tab');
type 'LOCATION\\DEB_IPM_sensitivity.txt'