syms phi(t) theta(t) gamma(t) h(t) N1(t) N2(t) Fat1u(t) Fat2u(t) Fat1s(t) Fat2s(t) r a b m mu1 mu2 g

psi = asin((sqrt(a^2-(r^2/cos(theta)^2)))/(a*sin(theta)));

rcm = [
    a*sin(psi)*sin(phi);
    -a*sin(psi)*cos(phi);
    h;
    ];

d1 = [
    a*cos(theta)*cos(phi + psi);
    a*sin(phi + psi);
    a*sin(theta)*cos(phi + psi);
    ];

d2 = - [
    a*cos(theta)*cos(phi - psi);
    a*sin(phi - psi);
    a*sin(theta)*cos(phi - psi);
    ];

p1 = simplify(rcm + d1);
p2 = simplify(rcm + d2);

n1v = [
    sin(psi)*sin(phi) + cos(theta)*cos(phi + psi);
    -sin(psi)*cos(phi) + sin(phi + psi);
    0;
    ];

n1 = n1v / norm(n1v);

n2v = [
    sin(psi)*sin(phi) - cos(theta)*cos(phi - psi);
    -sin(psi)*cos(phi) - sin(phi - psi);
    0;
    ];

n2 = n2v / norm(n2v);

zlinha = [
    -sin(theta);
    0;
    cos(theta);
    ];

t1 = cross(n1, zlinha);
t2 = cross(n2, zlinha);

Rz = [
    cos(phi), sin(phi), 0;
    -sin(phi), cos(phi), 0;
    0, 0, 1
    ];

Ry = [
    cos(theta), 0, sin(theta);
    0, 1, 0;
    -sin(theta), 0, cos(theta);
    ];

R = Rz * Ry;
Rinv = simplify(inv(R));

Ilinha = [
    m*(a^2 + b^2)/4, 0, 0;
    0, m*(a^2 + b^2)/4, 0;
    0, 0, m*(a^2 + b^2)/2;
    ];

D = [
    (a*sin(psi)*cos(phi))^2 + h^2, ((a*sin(psi))^2)*sin(phi)*cos(phi), -a*h*sin(psi)*sin(phi);
    ((a*sin(psi))^2)*sin(phi)*cos(phi), (a*sin(psi)*sin(phi))^2 + h^2, a*h*sin(psi)*cos(phi);
    -a*h*sin(psi)*sin(phi), a*h*sin(psi)*cos(phi), (a*sin(psi))^2;
    ];

I = simplify(simplify(R*Ilinha*Rinv) + D);

omega = [
    -diff(theta, t)*cos(theta)*sin(phi);
    diff(theta, t)*cos(theta);
    diff(phi, t) - diff(theta, t)*sin(theta)*sin(phi);
    ];

omegalinha = diff(gamma, t) * [
     -sin(theta);
    0;
    cos(theta);
    ];

W = [
    0;
    0;
    -m*g;
    ];

L = simplify(I*omega + Ilinha*omegalinha);
Tu = simplify(diff(L, t)); % unsteady
Ts = subs(Tu, [diff(theta, t), diff(theta, t, t), diff(h, t, t)], [0, 0, 0]); % steady

Tsys = cross(rcm, W) + N1 * cross(p1, n1) + N2 * cross(p2, n2) + mu1 * N1 *cross(p1, t1) + mu2*N2; % Falta o atrito (precisa determinar as relações de escorregamento)

P = m * diff(rcm, t);
Fu = simplify(diff(P, t)); % unsteady
Fs = simplify(subs(Fu, [diff(theta, t), diff(theta, t, t), diff(h, t, t)], [0, 0, 0])); % steady

Fsys = N1*n1 + N2*n2 + W + mu1*N1*t1 + mu2*N2*t2; % Falta o atrito (precisa determinar as relações de escorregamento)

eq1s = Ts == Tsys;
eq2s = Fs == Fsys;
% eq3s = Ilinha * [0; 0; diff(gamma, t);] == Coisa; % Precisa deteminar melhor, essa equação é os torques no frame linha

V1 = odeToVectorField(eq1s);
