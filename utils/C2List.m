function [List] = C2List(C)
NF = size(C,3);
List = [];
    for i=1:NF
        tp = C(1:NF,1:NF,i)*i;
        [x,y,z] = find(tp);
        List=vertcat(List,[x,y,z]);
    end

end