function out = fullContrast_subset(choice,P,CL,CR,uniqueC)
%ZL and ZR functions for the full contrast model

uniqueCL = unique(uniqueC(:,1));
uniqueCR = unique(uniqueC(:,2));
uniqueCL = uniqueCL(uniqueCL>0);
uniqueCR = uniqueCR(uniqueCR>0);

numPperZ = 1+length(uniqueCL)+length(uniqueCR);

numPZL = 1+length(uniqueCL);
numPZR = 1+length(uniqueCR);

switch(choice)
    case 'L'
        PZL = P(1:numPZL);
        
        N=length(CL);
        out = nan(N,1);
        for n=1:length(CL)
            out(n,1) = PZL*[1; CL(n)==uniqueCL];
        end
        
    case 'R'
        PZR = P(numPZL+1:(numPZL+numPZR));
        N=length(CL);
        out = nan(N,1);
        for n=1:length(CR)
            out(n,1) = PZR*[1; CR(n)==uniqueCR];
        end
        
    case 'paramLabels'
        ZLlabs = cell(numPZL,1);
        ZRlabs = cell(numPZR,1);
        
        ZLlabs{1} = 'Offset_L';
        ZRlabs{1} = 'Offset_R';
        
        ZLlabs(2:end) = cellstr(strcat(num2str(uniqueCL),'L_L'));
        ZRlabs(2:end) = cellstr(strcat(num2str(uniqueCR),'R_R'));
        
        out = [ZLlabs;ZRlabs]';
    case 'paramBounds'
        out=[-inf(1,numPZL+numPZR);inf(1,numPZL+numPZR)];
end

end