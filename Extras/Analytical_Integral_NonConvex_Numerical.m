function [ val ] = Analytical_Integral_NonConvex_Numerical( fnNumber, y, a , b )

val = 0;
switch fnNumber
    case 1
        %unit function
        f = @(x) y(x);
    case 2
        %x
        f = @(x) x.*y(x);
    case 3
        %y
        f = @(x) y(x).^2 / 2;
    case 4
        %x^2
        f = @(x) x.^2.*y(x);
    case 5 
        %xy
        f = @(x) (x.*y(x).^2)/2;
    case 6 
        %y^2
        f = @(x) (y(x).^3) / 3;
    case 7
        %x^3
        f = @(x) x.^3 .* y(x); 
    case 8
        %x^2 y
        f = @(x) (x.^2.*y(x).^2) / 2;
    case 9
        %xy^2
        f = @(x) (x.*y(x).^3) / 3;
    case 10
        %y^3
        f = @(x) (y(x).^4 / 4 );
    otherwise
        disp('other value')
end
        val = integral( f , a , b );
end