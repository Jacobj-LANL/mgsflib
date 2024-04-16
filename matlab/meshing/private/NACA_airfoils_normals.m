function NACA_airfoils_normals(digits, scale_pts, np)
if nargin < 3
    np = 50;
end
if nargin < 2
    scale_pts = [.2,.8,0.5]; % third value is y center
end
if nargin < 1
    digits = '0030';
end

m = str2num(digits(1))/100;
p = str2num(digits(2))/10;
t = str2num(digits(3:4))/100;

a0 = 0.2969;
a1 = 0.1260;
a2 = 0.3516;
a3 = 0.2843;
a4 = 0.1015;

syms x real
D = 5*t*(a0*sqrt(x) - a1*x - a2*x^2 + a3*x^3 - a4*x^4);
y = matlabFunction(D);
nrms = matlabFunction([-1,1/diff(D,x)]);

h = 1/(np-1);
t = [linspace(0,1,np)'; linspace(1-h,h,np-2)'];
nv = length(t);
normals = zeros(nv,2);
xs = zeros(nv,2);
for ii = 1:nv
    nn = nrms(t(ii));
    normals(ii,:) = nn/norm(nn);
    xs(ii,1) = t(ii);
    xs(ii,2) = y(t(ii));
    if ii > np
        xs(ii,2) = -xs(ii,2);
        normals(ii,2) = -normals(ii,2);
    end
    if t(ii) < 0.3
        normals(ii,:) = -normals(ii,:);
    end
    if t(ii) == 1.0
        normals(ii,:) = [-1,0];
    end
end
normals(1,1) = 1;

normals = -normals;
figure;
hold on
axis([-0.5,1.5,-0.5,1.5])
plot(xs(:,1),xs(:,2))
quiver(xs(:,1),xs(:,2),normals(:,1),normals(:,2))
end