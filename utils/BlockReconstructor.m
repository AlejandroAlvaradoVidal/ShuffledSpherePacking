function [out] = BlockReconstructor(I,NF,x,y)
cant=size(I,4);
z=size(I,3);
out=[];
c=1;
for i=1:x/NF
    temp=[];
    for j=1:y/NF
        block=I(:,:,:,c);
        temp=horzcat(temp,block);
        c=c+1;
    end
        out=vertcat(out,temp);
end


end