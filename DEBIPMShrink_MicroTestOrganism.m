function [E_Y, SSD, PopQs] = DEBIPMShrink_MicroTestOrganism(E_Ymin,E_Ymax,step,MatrixSize,Lm,Lp,lowval,rb)

% input: 
% E_Ymin: minimum feeding level; suggestion: 0.5
% E_Ymax: maximum feeding level; suggestion: 1
% step: step size between E_Ymin and E_Ymax: suggestion: 50
% MatrixSize: number of size bins of matrix approximation of DEBIPM. suggestion: 200

% to run code with above suggested parameter values: [E_Y, SSD, PopQs] = DEBIPMShrink_MicroTestOrganism(0.5,1,50,200);

% DEB parameters for a made-up microorganism

j = 200; % defines matrix size
%Lm = 19.433; % maximum length
%Lp = 16; % length at division 
%lowval= 0.5*Lp; % size at division (=0.5*Lp); limits of matrix approximation
hival = Lm; % to define limits of matrix approximation
kappa = 1; % fraction energy allocation to respiration (as opposed to reproduction)
E_Ystdev = 0.1; % sigma(Y)
%rb = -2.27E-03; % von Bertalanffy growth rate
mu = 0.01; % mortality rate

E_Y=linspace(E_Ymin,E_Ymax,step); % Feeding level
       
        for j=1:length(E_Y); j/length(E_Y);                     
                        
            [S, R, G, D, meshpts] = BigMatrixShrink_MicroOrg(MatrixSize,lowval,hival,Lp,Lm,E_Y(j),E_Ystdev,rb,kappa,mu);
            kernelDEB = G*S + D*R;
           % calculate lambda, ssd, and rv
            [W,d] = eig(kernelDEB); lambda1 = diag(d); imax = find(lambda1==max(lambda1)); 
            V=conj(inv(W)); lambda = lambda1(imax); % population growth rate
            w1=W(:,imax); v1 = real(V(imax,:))';
            ssd = w1/sum(w1); ssd = ssd'; % stable stage distribution

            % calculation of LRS
            [m,p]=size(kernelDEB); TS = G*S; Fmat = D*R; 
            R0mat = Fmat*inv((eye(m)-TS)); meanLRS = max(eig(R0mat));

            %  generation time
            GT = log(meanLRS)/log(lambda);
            
            if size(GT,2) ~= 1;
                GT = 0;
            end
            
            if size(lambda',2) ~= 1;
                lambda = 0;
            end
            
            PopQstemp = [E_Y(j),lambda',meanLRS,GT];
            PopQs(j,:) = PopQstemp;
            SSD(:,j) = ssd;
        end
        
        subplot(3,1,1); plot(E_Y,PopQs(:,2),'k-','LineWidth',2);  axis square; xline(0.9,'color','r'); yline(1.17,':'); text(0.87,0.535,"0.9",'color','r','FontSize',8); text(1.05,1.17,"1.17",'FontSize',8); box off
        ylabel('Population growth rate','FontSize',12); xlabel(' ');  set(gca,'FontSize',8)
        subplot(3,1,2); plot(E_Y,PopQs(:,3),'k-','LineWidth',2);   axis square; xline(0.9,'color','r'); yline(1.93,':'); text(0.87,-0.21,"0.9",'color','r','FontSize',8); text(1.05,1.93,"1.93",'FontSize',8); box off
        ylabel('R_0','FontSize',12);  xlabel(' ');  set(gca,'FontSize',8); 
        subplot(3,1,3); plot(E_Y,PopQs(:,4),'k-','LineWidth',2); axis square; xline(0.9,'color','r'); yline(4.24,':'); text(0.87,-6.1,"0.9",'color','r','FontSize',8);text(1.05,4.24,"4.24",'FontSize',8); box off;
        ylabel('Generation time','FontSize',12);  xlabel('Feeding level');  set(gca,'FontSize',8);

