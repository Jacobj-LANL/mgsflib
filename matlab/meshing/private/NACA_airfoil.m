function pts = NACA_airfoil(digits, scale_pts, np)

if nargin < 3
    np = 1000;
end
if nargin < 2
    scale_pts = [.2,.8,0.5]; % third value is y center
end
if nargin < 1
    digits = '0012';
end

m = str2num(digits(1))/100;
p = str2num(digits(2))/10;
t = str2num(digits(3:4))/100;

a0 = 0.2969;
a1 = 0.1260;
a2 = 0.3516;
a3 = 0.2843;
a4 = 0.1015;

h = 1/(np-1);
x = [linspace(0,1,np)'; linspace(1-h,h,np-2)'];
yt = @(x) 5*t*(a0*sqrt(x) - a1*x - a2*x^2 + a3*x^3 - a4*x^4);
yc = @(x) [(m/(p^2))*(2*p*x - x^2); (m/((1-p)^2))*((1-2*p) + 2*p*x - x^2)];
dyc_dx = @(x) [(2*m / (p^2))*(p-x); (m/((1-p)^2))*(p-x)];
nv = length(x);

% evaluating the functions
pts = zeros(nv,2);
if m == 0 && p == 0
    for ii = 1:np
        pts(ii,1) = x(ii);
        pts(ii,2) = yt(x(ii));
    end

    for ii = (np+1):nv
        pts(ii,1) = x(ii);
        pts(ii,2) = -yt(x(ii));
    end

else
    for ii = 1:np
        f = (x(ii)>p)+1;
        dd = dyc_dx(x(ii));
        theta = atan(dd(f));
        pts(ii,1) = x(ii) - yt(x(ii))*sin(theta);
        y = yc(x(ii));
        pts(ii,2) = y(f) + yt(x(ii))*cos(theta);
    end

    for ii = (np+1):nv
        f = (x(ii)>p)+1;
        dd = dyc_dx(x(ii));
        theta = atan(dd(f));
        pts(ii,1) = x(ii) + yt(x(ii))*sin(theta);
        y = yc(x(ii));
        pts(ii,2) = y(f) - yt(x(ii))*cos(theta);
    end
end

% scale the array
min_x = min(pts(:,1));
max_x = max(pts(:,1));

ratio = (max_x - min_x) / (scale_pts(2) - scale_pts(1));

pts = ((pts-min_x)/ratio) + [scale_pts(1),scale_pts(3)];
end