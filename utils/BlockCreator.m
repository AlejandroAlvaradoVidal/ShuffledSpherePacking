function [out] = BlockCreator(I,NF,NF2)
test={};
out=[];
c=1;
for i=1:size(I,1)/NF
    for j=1:size(I,2)/NF
        temp=I((i-1)*NF+1:i*NF,(j-1)*NF+1:j*NF,1:NF2);
        %out=cat(4,out,temp);
        test{c}=temp;
        c=c+1;
    end
    %out=cat(4,out,temp);
end

for i=1:c-1
    out=cat(4,out,test{i});
end
end