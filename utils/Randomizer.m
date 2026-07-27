function [C,info] = Randomizer(C,R)
NF=size(C,3);
mov=randi([-1,1],NF,NF);
Cold=C;
C=zeros(size(C));
for i=1:NF
    for j=1:NF
        temp=squeeze(Cold(i,:,j));
        temp=circshift(temp,mov(i,j));
        C(i,:,j)=temp;
    end
end
%%
pat=ones(size(C));
for i=1:NF
    pat(:,:,i)=i.*pat(:,:,i);
end
x=C.*pat;
d=ComputedistanceGen(x,NF);%min([d1,d2,d3]);

sphere = 0;
for i =1:nnz(C)
    sphere = sphere + (4/3).*pi.*((d/2).^3);
end
cube = (NF+1).^3;
density= gather(sphere/cube);
%disp("Sphere Packing Density "+ density + " diameter "+ d);
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
info=[Trans,Cond,diameter,density];
end