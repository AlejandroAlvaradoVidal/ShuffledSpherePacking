function [C] = List2C(List)
NF = max(List(:,3));
C=zeros(NF,NF,NF);
    for i=1:size(List,1)
        id=List(i,:);
        x=id(1);
        y=id(2);
        z=id(3);
        C(x,y,z)=1;
    end
end