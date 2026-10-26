function [i1,i2] = voigt2xyz(I)

switch I
    case 1
        i1 = 1; i2 = 1;
    case 2
        i1 = 2; i2 = 2;
    case 3
        i1 = 3; i2 = 3;
    case 4
        i1 = 2; i2 = 3;
    case 5
        i1 = 1; i2 = 3;
    case 6
        i1 = 1; i2 = 2;
end

end