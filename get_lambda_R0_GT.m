function [lambda, R0, GT] = get_lambda_R0_GT(All_Species,E_Y,i)
   
    
    index = find(round(All_Species(i,:,1),2) == E_Y);
    lambda = All_Species(i,index,2);
    R0 = All_Species(i,index,3);
    GT = All_Species(i,index,4);