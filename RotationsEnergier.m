% Variabler

rot_energinivaer = 75;
vib_energinivaer = 81;

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
D_e = 6.12147e-6*100;
%omega = 2169.81358*2*pi*2.99792458*10^(10);
%omega_xe = 13.28831*2*pi*2.99792458*10^(10);
omega = 2169.81358*100;
omega_xe = 13.28831*100;
% Konstanter

hbar = 1.0546e-34;
h = 6.6261e-34;
c = 299792458;

M = M_A*M_B/(M_A + M_B);
D = 2*pi*hbar*c*D_e;
k = omega^2*M;

a = sqrt(2*M*(omega_xe/(2*pi*c))/hbar);
%D_w = (M*omega^2)/(a^2*2*h*c);
D_w = 1.7800979e-18;
a = (omega*2*pi*c)/sqrt(2*D_w/M);

% Rotation

fprintf('Rotationsenergier:\nn\tE_n\t\tlambda\n');

E_rot_korr      = zeros(1,rot_energinivaer);
E_rot_korr_distance = zeros(1,rot_energinivaer-1);
%E_rot_no_korr   = zeros(1,rot_energinivaer);

for (n = 0:rot_energinivaer)
    E_rot_korr(n+1)          = (n*(n+1)*hbar^2)/(2*M*R_e^2) - D*n^2*(n+1)^2;
    
    if n > 0
        E_rot_korr_distance(n) = E_rot_korr(n+1) - E_rot_korr(n);
    end
    
    %E_rot_no_korr(n+1)  = (n*(n+1)*hbar^2)/(2*M*R_e^2);
    lambda = h*c/E_rot_korr(n+1);
    
    %fprintf('%i\t%e\t%e\n', n, E_rot_korr(n+1), lambda);
end 

fprintf('\nMedelavstand rot eng: %e    Vaglangd: %e\n', mean(E_rot_korr_distance)./(1.6022e-19), h*c/mean(E_rot_korr_distance));
fprintf('Maxavstand rot eng:   %e    Vaglangd: %e\n', max(E_rot_korr_distance)./(1.6022e-19), h*c/max(E_rot_korr_distance));

% Vibration (morse)

fprintf('Vibrationsenergier\nn\tE_n\t\tlambda\n');

E_vib = zeros(1,vib_energinivaer);
E_vib_distance = zeros(1, vib_energinivaer-1);

for (n = 0:vib_energinivaer)
    %E_vib(n + 1) = hbar*((0.5 + n)*omega - (0.5 + n)^2*omega_xe);
    E_vib(n + 1) = h*c*((0.5 + n)*omega - (0.5 + n)^2*omega_xe);
    lambda = h*c/(E_vib(n+1) - E_vib(1));
    
    if n > 0
        E_vib_distance(n) = E_vib(n+1) - E_vib(n);
    end
    
    fprintf('%i\t%e\t%e\n', n, E_vib(n+1), lambda); 
end

n_max = (2*D_w - h*c*omega)/(h*c*omega);

fprintf('n_max: %e', n_max);
fprintf('\nMedelavstand vib eng: %e    Vaglangd: %e\n', mean(E_vib_distance)./(1.6022e-19), h*c/mean(E_vib_distance));
fprintf('Maxavstand vib eng:   %e    Vaglangd: %e\n', max(E_vib_distance)./(1.6022e-19), h*c/max(E_vib_distance));
fprintf('Minavstand vib eng:   %e    Vaglangd: %e\n', min(E_vib_distance)./(1.6022e-19), h*c/min(E_vib_distance));

% Plot all the things

% Morse VS harmonisk
% figure;
x = linspace(-5e-10,5e-10, 1000);
%y_morse = h*c*D_w*(1-exp(-a*x)).^2;
y_morse = D_w*(1-exp(-a*x)).^2;
y_harm = 0.5*k*x.^2;
P = InterX([x;y_morse], [x;E_vib(end)*ones(1,length(x))]);
% 
% plot(x/(1e-10), y_morse/(1.6022e-19), 'b', 'LineWidth', 1);
% hold on;
% plot(x/(1e-10), y_harm/(1.6022e-19), 'r', 'LineWidth', 1);
% 
% box on;
% grid on;
% ylim([0 2]);
% xlabel('$R - R_e$ (\AA)', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('Potentiell energi (eV)', 'Interpreter', 'latex', 'FontSize', 18);
% title('Harmoniska- och Morse-potentialen vid tv{\aa}atomig vibration', 'Interpreter', 'latex', 'FontSize', 18);
% h = legend('Morse', 'Harmonisk', 'Location', 'SouthEast');
% set(h, 'Interpreter', 'latex');
% h.FontSize = 18;

%  Potential
figure;
p_pot = plot(x/(1e-10), y_morse/(1.6022e-19), 'k', 'LineWidth', 2);
hold on;

%  Rot
for(n=0:rot_energinivaer)
    P = InterX([x;y_morse], [x;E_rot_korr(n+1)*ones(1,length(x))]);
    if ~isempty(P)
        x1 = linspace(P(1,1),P(1,2))./(1e-10);
        %x1 = linspace(-0.01,0.01, 100);
        y1 = E_rot_korr(n+1)*ones(1,length(x1));
        p_rot = plot(x1, y1/(1.6022e-19), 'b');
    end
end

%  Vib
for (n=0:vib_energinivaer)
    P = InterX([x;y_morse], [x;E_vib(n+1)*ones(1,length(x))]);
    x1 = linspace(P(1,1),P(1,2), 1000);
    y1 = ones(1,length(x1))*E_vib(n+1);
    p_vib = plot(x1/(1e-10),y1/(1.6022e-19), 'r', 'LineWidth', 2);
end

% Plot settings
ylim([0 (E_vib(end) + 3e-20)/(1.6022e-19)]);
xlim([(P(1,1)-0.05e-10)/(1e-10) (P(1,2)+0.05e-10)/(1e-10)]);

xlabel('$R - R_e$ (\AA)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Energi (eV)', 'Interpreter', 'latex', 'FontSize', 18);

h = legend([p_rot p_vib p_pot], {'$E_{rot}(n)$', '$E_{vib}(n)$', '$V_{morse}(R)$'}, 'Location', 'SouthEast');
set(h, 'Interpreter', 'latex');
h.FontSize = 18;

box on;
%grid on;

title('Vibrations- och rotationsenerginiv{\aa}er f\"or CO', 'Interpreter', 'latex', 'FontSize', 18);

% Only vib
figure;

p_pot = plot(x/(1e-10), y_morse/(1.6022e-19), 'k', 'LineWidth', 2);
hold on;
for (n=0:vib_energinivaer)
    P = InterX([x;y_morse], [x;E_vib(n+1)*ones(1,length(x))]);
    x1 = linspace(P(1,1),P(1,2), 1000);
    y1 = ones(1,length(x1))*E_vib(n+1);
    p_vib = plot(x1/(1e-10),y1/(1.6022e-19), 'r', 'LineWidth', 1);
end

% Plot settings
ylim([0 15]);
xlim([(P(1,1)-0.25e-10)/(1e-10) (P(1,2)+0.25e-10)/(1e-10)]);

xlabel('$R - R_e$ (\AA)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Energi (eV)', 'Interpreter', 'latex', 'FontSize', 18);

h = legend([p_vib p_pot], {'$E_{vib}(n)$', '$V_{morse}(R)$'}, 'Location', 'SouthEast');
set(h, 'Interpreter', 'latex');
h.FontSize = 18;

box on;
%grid on;

title('Vibrationsenerginiv{\aa}er f\"or CO', 'Interpreter', 'latex', 'FontSize', 18);


% Plot korr och icke-korr f?r rot

% figure;
% 
% x_korr      = linspace(-1e-10,-0.1e-10, 100) ./ (1e-10);
% 
% for(n=0:rot_energinivaer)
%     y_korr = E_rot_korr(n+1)*ones(1,length(x_korr));
%     p_rot = plot(x_korr, y_korr/(1.6022e-19), 'b');
%     hold on;
% end
% 
% %grid on;
% box on;
% 
% xlabel('Enhetsl\"os', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('Energi (eV)', 'Interpreter', 'latex', 'FontSize', 18);
% title('Centrifugalkorrigerade rotationsenerginiv{\aa}er f\"or CO', 'Interpreter', 'latex', 'FontSize', 18);

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


