% Variabler

rot_energinivaer = 20????;
vib_energinivaer = 20;

% H2
%M_A = 1.6737236e-27;
%M_B = 1.6737236e-27;
%R_e = 0.74144e-10;
%D_e = 0.0471*1e2;
%omega = 4401.21*2.99792458*10^(10);
%omega_xe = 121.33*2.99792458*10^(10);

% CO
M_A = 1.9944235e-26;
M_B = 2.6566962e-26;
R_e = 1.128323e-10;
D_e = 6.12147e-6*1e2;
omega = 1743.41*2.99792458*10^(10);
omega_xe = 14.36*2.99792458*10^(10);

% Konstanter

hbar = 1.0546e-34;
h = 6.6261e-34;
c = 299792458;

M = M_A*M_B/(M_A + M_B);
D = 2*pi*hbar*c*D_e;
k = omega^2*M;

a = sqrt(2*M*omega_xe/hbar);
D_w = (M*omega^2)/(a^2*2*h*c);

% Rotation

fprintf('Rotationsenergier:\nn\tE_n\t\tlambda\n');

E_rot = zeros(1,rot_energinivaer);

for (n = 0:rot_energinivaer)
    E_rot(n+1) = (n*(n+1)*hbar^2)/(2*M*R_e^2) - D*n^2*(n+1)^2;
    lambda = h*c/E_rot(n+1);
    
    fprintf('%i\t%e\t%e\n', n, E_rot(n+1), lambda);
end

% Vibration (morse)

fprintf('Vibrationsenergier\nn\tE_n\t\tlambda\n');

E_vib = zeros(1,vib_energinivaer);

for (n = 0:vib_energinivaer)
    E_vib(n + 1) = (0.5 + n)*hbar*omega - (0.5 + n)^2*hbar*omega_xe;
    lambda = h*c/(E_vib(n+1) - E_vib(1));
    
    fprintf('%i\t%e\t%e\n', n, E_vib(n+1), lambda); 
end

% Plot all the things

%  Potential
x = linspace(-5e-10,5e-10, 1000);
y = h*c*D_w*(1-exp(-a*x)).^2;
P = InterX([x;y], [x;E_vib(end)*ones(1,length(x))]);

p_pot = plot(x/(1e-10), y/(1.6022e-19), 'k', 'LineWidth', 2);
hold on;

%  Rot
for(n=0:rot_energinivaer)
    P = InterX([x;y], [x;E_rot(n+1)*ones(1,length(x))]);
    if ~isempty(P)
        x1 = linspace(P(1,1),P(1,2))./(1e-10);
        y1 = E_rot(n+1)*ones(1,length(x1));
        p_rot = plot(x1, y1/(1.6022e-19), 'b');
    end
end

%  Vib
for (n=0:vib_energinivaer)
    P = InterX([x;y], [x;E_vib(n+1)*ones(1,length(x))]);
    x1 = linspace(P(1,1),P(1,2), 1000);
    y1 = ones(1,length(x1))*E_vib(n+1);
    p_vib = plot(x1/(1e-10),y1/(1.6022e-19), 'r', 'LineWidth', 2);
end

% Plot settings
ylim([0 (E(end) + 1e-20)/(1.6022e-19)]);
xlim([(P(1,1)-0.25e-10)/(1e-10) (P(1,2)+0.25e-10)/(1e-10)]);

xlabel('$\Delta R$ (\AA)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Energi (eV)', 'Interpreter', 'latex', 'FontSize', 18);

h = legend([p_rot p_vib p_pot], {'$E_{rot}(n)$', '$E_{vib}(n)$', 'Potential $V(\Delta R)$'}, 'Location', 'SouthEast');
set(h, 'Interpreter', 'latex');
h.FontSize = 18;

box on;
grid on;

title('Vibrations- och rotationsenerginiv{\aa}er', 'Interpreter', 'latex', 'FontSize', 18);

%% Vibration (kvadratisk approx.)

energinivaer = 10;

M_A = 1.6737236e-27;
M_B = 1.6737236e-27;
%omega = 4161.166*2.99792458*10^(10);
%omega = 2321.4*2.99792458*10^(10);
omega = 4401.21*2.99792458*10^(10);

hbar = 1.0546e-34;
h = 6.6261e-34;
c = 299792458;

M = M_A*M_B/(M_A + M_B);
k = omega^2*M;

fprintf('Vibrationsenergier\nn\tE_n\t\tlambda\n');

E = zeros(1,energinivaer);

for (n = 0:energinivaer)
    E(n + 1) = (0.5 + n)*hbar*omega;
    lambda = h*c/E(n+1);
    
    fprintf('%i\t%e\t%e\n', n, E(n+1), lambda); 
end

x = [-5e-10:1e-12:5e-10];
y = 0.5*k*x.^2;

[diff, index] = min(abs(y - E(end)));
x_width = abs(x(index));

plot(x, y);
hold on;
ylim([0 E(end) + 1e-20]);
xlim([-x_width-0.25e-10 x_width+0.25e-10]);


for (n=0:energinivaer)
    [diff, index] = min(abs(y - E(n+1)));
    x_width = abs(x(index));
    x1 = [-x_width:1e-12:x_width];
    y1 = ones(1,length(x1))*E(n+1);
    plot(x1,y1);
end


