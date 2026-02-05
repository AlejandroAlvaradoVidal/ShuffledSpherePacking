function [C] = Cat_OSP2(NF,M,R)
NF2=ceil(NF/M);
NFo=NF;
if NF2*M<NF
    NF2=ceil(NF/M)+1;
    NF=M*NF2;
end
[G] = RunEQ10MUX2(NF2,2*M);
id=round(linspace(1,2*M,M)); 
G=G(:,:,id);
C=[];
for i=1:M
    C=cat(3,C,G2C(G(:,:,i)));
end
C=repmat(C,[M,M]);
C=C(1:NFo,1:NFo,1:NFo);
Cond=GetCond(C);
Trans=sum(C,"all")/NF^3;
disp("Transmittance "+ Trans + " Condition Number "+ Cond);
%disp("Mux "+ num2str(unique(sum(C,3))));
%%
if R==1
    [density,diameter]= ComputeDensityIrregularSP(C);
    disp("Irregural: Density "+ density + " Diameter "+ diameter);
else
    [density,diameter]= ComputeDensityRegularSP(C);
    disp("Regural: Density "+ density + " Diameter "+ diameter);
end
end