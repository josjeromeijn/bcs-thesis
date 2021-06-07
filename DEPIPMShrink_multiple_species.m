% This code is to run the DEB-IPM over multiple species. 8 parameters are
% needed:

warning('off')

% E_Ymin: minimum feeding level; suggestion: 0.5
% E_Ymax: maximum feeding level; suggestion: 1
% step: step size between E_Ymin and E_Ymax: suggestion: 50
% MatrixSize: number of size bins of matrix approximation of DEBIPM. suggestion: 200
% Lm,Lp,lowval (Lb) and rb are retrieved from local .txt files 

E_Ymin = 0.5;
E_Ymax = 1;
step = 51;
MatrixSize = 200;

%Read in all Lm of species, adjust location to path of files
lm_table = readtable('LOCATION\\lm.txt'); 

%Read in all vB growth rates of species, adjust location to path of files
vB_table = readtable('LOCATION\\growthrates.txt');

%Create list of all species id's
all_id = lm_table.id;

%Create loop over all species id's
for i = [1:length(all_id)]
    id = all_id(i);
    
    %Find Lm of id
    Lm = lm_table.Lm(i);
    
    %Find file of id, adjust path to location of files
    file_name = 'LOCATION\\VonBertalanffy\\id' + string(id) + '.txt';
    single_id_table = readtable(file_name);
    
    %Find lowval (Lb) and Lp of id
    lowval = single_id_table.size_t(1);
    Lp = single_id_table.size_t(end);
    
    %Find rb of id
    rb = vB_table.Var2(i);
    
    %Calculate SSD & PopQs for each value of E_Y
    [E_Y, SSD, PopQs] = DEBIPMShrink_MicroTestOrganism(E_Ymin,E_Ymax,step,MatrixSize,Lm,Lp,lowval,rb);
        
    %Save plots of pop growth rate, R0 and generation time of id
    % file_name = 'LOCATION\\id' + string(id) + '_E_Ymin_' + string(E_Ymin) + '.jpg';
    % saveas(gcf,file_name)
    
    %Store PopQs (50x4 = E_Y, pop growth rate, R0, generation time). This
    %matrix contains all pop growth rates, R0 and generation times for each
    %value of E_Y, per species id. (=17x51x4)
    All_Species_0_5(i,:,:) = PopQs;
    
end
