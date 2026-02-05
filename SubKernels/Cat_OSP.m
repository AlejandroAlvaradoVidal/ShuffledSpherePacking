function [C] = Cat_OSP(NF,M,R)
[G] = RunEQ10MUX(ceil(NF/M),2*M);
id=round(linspace(1,2*M,M)); 
G=G(:,:,id);
C=[];
for i=1:M
    C=cat(3,C,G2C(G(:,:,i)));
end
C=repmat(C,[M,M]);
C=C(1:NF,1:NF,1:NF);
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