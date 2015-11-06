function out = fullContrast(choice,P,CL,CR,uniqueC)
%ZL and ZR functions for the full contrast model

uniqueCL = unique(uniqueC(:,1));
uniqueCR = unique(uniqueC(:,2));
uniqueCL = uniqueCL(uniqueCL>0);
uniqueCR = uniqueCR(uniqueCR>0);

numPperZ = 1+length(uniqueCL)+length(uniqueCR);

switch(choice)
    case 'L'
        PZL = P(1:numPperZ);
        
        N=length(CL);
        out = nan(N,1);
        for n=1:length(CL)
            out(n,1) = PZL*[1; CL(n)==uniqueCL; CR(n)==uniqueCR];
        end
    case 'R'
        PZR = P(numPperZ+1:2*numPperZ);
        N=length(CL);
        out = nan(N,1);
        for n=1:length(CR)
            out(n,1) = PZR*[1; CL(n)==uniqueCL; CR(n)==uniqueCR];
        end
    case 'paramLabels'
        out = cell(numPperZ,2);
        out(1,:) = {'Offset_L','Offset_R'};
        out(2:end,1) = cellstr(strcat([strcat(num2str(uniqueCL),'L'); strcat(num2str(uniqueCR),'R')],'_L'));
        out(2:end,2) = cellstr(strcat([strcat(num2str(uniqueCL),'L'); strcat(num2str(uniqueCR),'R')],'_R'));
        
        out = out(:)';
    case 'paramBounds'
        out=[-inf(1,2*numPperZ);inf(1,2*numPperZ)];
end

end