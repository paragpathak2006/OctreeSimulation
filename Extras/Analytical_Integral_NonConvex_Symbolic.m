function [ val ] = Analytical_Integral_NonConvex_Symbolic( fnNumber, y, a , b )

syms x

val = 0;
switch fnNumber
    case 1
        %unit function
        val = int( y,a,b);
    case 2
        %x
        val = int( x*y, a,b );
    case 3
        %y
        val = int( y^2 / 2 , a , b );
    case 4
        %x^2
        val = int( x^2.*y, a, b);
    case 5 
        %xy
        val = int( x*y^2/2, a,b);
    case 6 
        %y^2
        val = int( y^3 / 3, a, b);
    case 7
        %x^3
        val = int( x^3 * y , a , b );
    case 8
        %x^2 y
        val = int( x^2*y^2 / 2, a , b );
    case 9
        %xy^2
        val = int( x*y^3 / 3, a , b );
    case 10
        %y^3
        val = int( y^4 / 4, a , b );
    otherwise
        disp('other value')
end

end