rv = 1.5;
Rv = 2;
ru = 1.6;
Ru = 3;
%
u0 = sqrt(Ru^2-ru^2);
theta0 = acos(ru/Ru);
d0 = Rv-rv;
%
L = sqrt(d0^2/2)*1.8;
alpha10 = acos(d0/(2*L));
%
u = linspace(2.25,3.22);
%
theta = acos((ru^2 + Ru^2 - u.^2)/2/ru/Ru);
dtheta = theta-theta0;
d = sqrt(rv^2 + Rv^2-2*rv*Rv*cos(dtheta));
psi = acos((d.^2 + Rv^2 - rv^2)./(2*d*Rv)).*sign(dtheta);
phi = acos(d/2/L);
%
alpha1 = phi - psi;
%
plot(u,alpha1*180/pi);
xlabel('$u$ [m]','interpreter','latex')
ylabel('$\alpha_1$ [${}^\circ$]','interpreter','latex')
title('Guide vane angle $\alpha_1$ as a function of actuator position $u$',...
    'interpreter','latex')
grid on
hold on
plot(u0,alpha10*180/pi,'ko')
hold off