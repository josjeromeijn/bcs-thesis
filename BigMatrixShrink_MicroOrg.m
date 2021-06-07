
%%%%%%%%%%%%%%%% The 'big matrix' M of size n x n %%%%%%%%%%%%%%%%%%%%%%%%%
function [S, R, G, D, y] = BigMatrixShrink_MicroOrg(n,lowval,hival,Lp,Lm,E_Y,E_Ystdev,rb,kappa,mu)
	
% upper and lower integration limits
	L = 0.9*lowval; U = 1.1*hival;
	% boundary points b and mesh points y
	c = 0:1:n; b = L+c.*(U-L)./n;
    y = 0.5.*(b(1:n)+b(2:n+1));
    
    % define probability of reproducing (step function): 
    Prob_reproduction = y>=Lp;
    
	% create S, R, G and D matrices
	S = SurvivalFunction(y,mu,Prob_reproduction,E_Y,kappa,Lm); S=repmat(S,n,1); S=S.*eye(n); 
    G = GrowthFunction(y,E_Y,E_Ystdev,Lm,rb); 
    R = ReproductionFunction(y,Prob_reproduction); R=repmat(R,n,1); R=R.*eye(n); 
	D = DevelopmentFunction(y,lowval); 
	% scale D and G so columns sum to 1
	[m,n]=size(G); G = G./repmat(sum(G),m,1);
    [m,n]=size(D); D = D./repmat(sum(D),m,1); 

% Compute the kernel component functions from the fitted models

% survival
function [sx] = SurvivalFunction(x,mu,Prob_reproduction,E_Y,kappa,Lm)
     sx1 =x<=Lm*E_Y/kappa; sx1=sx1*exp(-mu); % probability of survival of those that do not starve
     sx=sx1.*(1 - Prob_reproduction); % only individuals that do not divide, survive
    
% growth
function gxy = GrowthFunction(x,E_Y,E_Ystdev,Lm,rb) 
    for i=1:length(x)
        mux = x(i).*exp(-rb) + (1 - exp(-rb))*Lm*E_Y;
        sigmax = (1 - exp(-rb))*Lm*E_Ystdev;
        fac1 = sqrt(2.*pi).*sigmax;
        fac2 = ((x'-mux).^2)./(2.*sigmax.^2);
        gxy(:,i) = (exp(-fac2)./fac1);
    end

% reproduction
function [rx] = ReproductionFunction(x,Prob_reproduction)
    rx1 = ones(1,length(x))*2; 
    rx=rx1.*Prob_reproduction;% vector where reproducing inds divide into two offspring
    
% % development
% function [dxy] = DevelopmentFunction(x)
%         kidsize.mean = 0.5*x; % take half of mean size of individuals (as they divide)
%         kidsize.var2 = 0.0001; 
%         kidsize.var = sqrt(kidsize.var2);
%         fac1 = sqrt(2*pi).*kidsize.var;
%         fac2 = ((x'-kidsize.mean').^2)./(2.*kidsize.var2);
%         dxy = exp(-fac2)./fac1;
%     dxy=repmat(dxy,1,length(dxy));
    

function [dxy] = DevelopmentFunction(x,Lb)
    xtemp = x<=Lb; 
    dxy=repmat(xtemp',1,length(xtemp));
    dxy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%