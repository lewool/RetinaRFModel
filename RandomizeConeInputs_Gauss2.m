function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=RandomizeConeInputs_Gauss1(Em)

%Determine some parameters for centers and surrounds
CSRadRatio=6; %Let's say that a surround radius is ~6x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=.75; %Let's say that surround gain is ~3/4 that of the center (Croner & Kaplan, Vis. Res. 1995)

%How many cones comprise an RF center at any eccentricity?
ConesToCenter=ceil(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B

%Multiply that value by the C:S area ratio (if Rs=CSRadRatio*Rc, then
%As=(CSRadRatio^2)*Ac)
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

%Determine a 1-mm^2 patch of cones for a particular retinal eccentricity
ConesPerMM=ceil(19890*Em^(-0.6331)); %Curcio 1990

% Assign a retinal L:M:S ratio from which any one cell comes
% (L:M:S ratios vary across retinas, lognormally)
LogLMDistribution=makedist('Normal','mu',0.47,'sigma',0.74); %Mu and sigma computed from a fit of Dacey et al. (2000) JOSA 17:589-596, Fig 5A
RandLM=exp(random(LogLMDistribution));
pL=RandLM/(RandLM+1); %Likelihood of L
pM=1-pL; %Likelihood of M
pS=0; %Likelihood of S

%Assign the patch (a vector of length ConesPerMM) to L, M, and S
%'zones' (for assigning discrete cones later)
LLIndex=0;
LUIndex=ceil(ConesPerMM*pL);
MLIndex=LUIndex+1;
MUIndex=MLIndex+ceil(ConesPerMM*pM);
SLIndex=MUIndex+1;
SUIndex=ConesPerMM;

%Determine the dimensions of the NxM square matrix needed to place all
%the cells in a 2D arrangement
ConesToCenterDim=ceil(sqrt(ConesToCenter));
ConesToSurroundDim=ceil(sqrt(ConesToSurround));

%Compute a responsivity Gaussian for one dimension of the 2D grid

%Choose a value for the sigma (cone contribution falloff)
ApertureSigma=.34; %%%%Can't remember why we chose this value...

%Symmetric Gaussian correction
if rem(ConesToCenter,2)==1
    CenterStrength=normpdf(1:ConesToCenterDim,ceil(ConesToCenterDim/2),ApertureSigma*ConesToCenterDim);
else
    CenterStrength=normpdf(1:ConesToCenterDim,(ConesToCenterDim+1)/2,ApertureSigma*ConesToCenterDim);
end

if rem(ConesToSurround,2)==1
    SurroundStrength=normpdf(1:ConesToSurroundDim,ceil(ConesToSurroundDim/2),ApertureSigma*ConesToSurroundDim);
else
    SurroundStrength=normpdf(1:ConesToSurroundDim,(ConesToSurroundDim+1)/2,ApertureSigma*ConesToSurroundDim);
end

%Compute the outer product of 2 1D Gaussians to create a 2D Gaussian of the
%receptive field "center" and "surround" responsivity (there are no cones placed in the
%patch yet...)
CenterSurface=CenterStrength'*CenterStrength;
SurroundSurface=SurroundStrength'*SurroundStrength;

CenterWeights=sort(reshape(CenterSurface,[1,ConesToCenterDim^2]),'descend');
SurroundWeights=sort(reshape(SurroundSurface,[1,ConesToSurroundDim^2]),'descend');

%Now randomly assign L, M, and S cones to the surround
SurroundAssign=randi(ConesPerMM,1,ConesToSurround);

%Select a subsection of surround cones to serve as the center
%(the retina double dips)
CenterAssignIndex0=1;
CenterAssignIndexN=ConesToCenter;
CenterAssign=SurroundAssign(CenterAssignIndex0:CenterAssignIndexN);

%Determine which center values are Ls, Ms, and Ss...
CenterAssignL=CenterAssign<=LUIndex;
CenterAssignM=(CenterAssign>MLIndex & CenterAssign<=MUIndex);
CenterAssignS=(CenterAssign>SLIndex & CenterAssign<=SUIndex);

%Weight the cone contributions given the responsivity Gaussian computed
%above
cpL=(sum(CenterAssignL.*CenterWeights(1:ConesToCenter)))/(sum(CenterWeights(1:ConesToCenter)));
cpM=(sum(CenterAssignM.*CenterWeights(1:ConesToCenter)))/(sum(CenterWeights(1:ConesToCenter)));
cpS=(sum(CenterAssignS.*CenterWeights(1:ConesToCenter)))/(sum(CenterWeights(1:ConesToCenter)));

%Just keep track of the actual cone count
cnL=sum(CenterAssignL);
cnM=sum(CenterAssignM);
cnS=sum(CenterAssignS);

%Determine which surround values are Ls, Ms, and Ss...
SurroundAssignL=SurroundAssign<=LUIndex;
SurroundAssignM=(SurroundAssign>MLIndex & SurroundAssign<=MUIndex);
SurroundAssignS=(SurroundAssign>SLIndex & SurroundAssign<=SUIndex);

%Weight the cone contributions given the responsivity Gaussian computed
%above
spL=(sum(SurroundAssignL.*SurroundWeights(1:ConesToSurround)))/(sum(SurroundWeights(1:ConesToSurround)));
spM=(sum(SurroundAssignM.*SurroundWeights(1:ConesToSurround)))/(sum(SurroundWeights(1:ConesToSurround)));
spS=(sum(SurroundAssignS.*SurroundWeights(1:ConesToSurround)))/(sum(SurroundWeights(1:ConesToSurround)));

%Just keep track of the actual cone count
snL=sum(SurroundAssignL);
snM=sum(SurroundAssignM);
snS=sum(SurroundAssignS);

cnAll=[cnL cnM cnS];
snAll=[snL snM snS];

cpAll=[cpL cpM cpS];
spAll=[spL spM spS];