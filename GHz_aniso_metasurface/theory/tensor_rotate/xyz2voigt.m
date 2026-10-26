function I = xyz2voigt(i1,i2)
if     [i1,i2] == [1,1]
        I = 1;
elseif [i1,i2] == [2,2]
        I = 2;
elseif [i1,i2] == [3,3]
        I = 3;
elseif [i1,i2] == [2,3]
        I = 4;
elseif [i1,i2] == [3,2]
        I = 4;
elseif [i1,i2] == [1,3]
        I = 5;
elseif [i1,i2] == [3,1]
        I = 5;
elseif [i1,i2] == [1,2]
        I = 6;
elseif [i1,i2] == [2,1]
        I = 6;
end
end