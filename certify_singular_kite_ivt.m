function certify_singular_kite_ivt()
%CERTIFY_SINGULAR_KITE_IVT Certificates for the singular family and mass surface.
%
% Requires Symbolic Math Toolbox. Run from MATLAB with:
%     certify_singular_kite_ivt
%
% This is a single function file with local helper functions below.
% Symbolic integers perform exact arithmetic. Rational input endpoints are
% represented exactly; fixed-point interval products and reciprocals are
% rounded outward at 60 decimal places using symbolic floor/ceil.
% No floating-point tolerance or vpa value is used in a certificate check.
%
% Proof structure.
%   (1) An 18-by-18 exact rational subdivision certifies the geometric
%       domain, E_x<0, E_y<0, and the coefficient bounds on
%       [a,b] x [c,d] x {z_*}.
%   (2) The first IVT and E_y<0 give a unique continuous graph y=h(x)
%       satisfying E(x,h(x),z_*)=0.
%   (3) Certified signs of K3 at x=a and x=b give a common zero of E and
%       K3 by a second IVT.
%   (4) Uniform coefficient bounds and separated upper bounds certify the
%       full positive-mass interval I(p_*)=(0,tau(p_*)), where
%       tau=W/((alpha-1)*Delta6) and 0.0155<tau<0.0157.  In particular,
%       every 0<t<=1/100 is allowed.  The two endpoints have m3=0 and
%       m1=0, respectively, with the other four masses strictly positive.
%   (5) A nonzero (x,y)-Jacobian minor at the root gives the local analytic
%       curve by the implicit function theorem.  Extension of the strict
%       coefficient bounds to nearby z is by continuity after shrinking
%       the analytic neighborhood; it is not claimed to be a certified
%       three-dimensional interval box.
%   (6) Two independent dual-number passes give the x,y,z partials.
%       Implicit differentiation and the chain rule certify
%       179<partial_z S1(z_*,1/100)<197, where
%       S(z,t)=(m1(sigma(z),t),t,m5(sigma(z),t)).  Since m3=t, the
%       (m1,m3) projection has a nonzero (z,t)-Jacobian determinant.
%       The analytic inverse function theorem gives a local embedded mass
%       surface; no global embedding or injectivity is asserted.
%
% Displayed decimals are rounded only for readability. All assertions use
% the stored 60-place outward interval endpoints as exact symbolic integers.
if isempty(ver('symbolic'))
    error('certify_singular_kite_ivt:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
C.prec = 60;
C.subdivisions = 18;
C.scale = sym(10)^C.prec;
C.scale2 = C.scale^2;
Z = rat('2783/5000');
XLEFT = rat('0.9780199');
XRIGHT = rat('0.9780201');
YLOW = rat('0.2957641');
YHIGH = rat('0.2957644');
% Exact rational enclosure of sqrt(3).  The square comparisons below are
% exact and prove SQRT3LOW < sqrt(3) < SQRT3HIGH.
SQRT3LOW = rat('173205/100000');
SQRT3HIGH = rat('173206/100000');
require(ratLT(ratMul(SQRT3LOW, SQRT3LOW), rat('3')));
require(ratLT(rat('3'), ratMul(SQRT3HIGH, SQRT3HIGH)));
sqrt3 = iv(SQRT3LOW, SQRT3HIGH, C);
% Endpoint enclosures for the two zeros of E(x,.,z_*).  These are the
% intervals printed in Appendix D of the manuscript.
brackets(1).xp = XLEFT;
brackets(1).yl = rat('0.2957643031163');
brackets(1).yr = rat('0.2957643031164');
brackets(1).kl = rat('7/10000000000');
brackets(1).kr = rat('8/10000000000');
brackets(2).xp = XRIGHT;
brackets(2).yl = rat('0.29576420207515');
brackets(2).yr = rat('0.29576420207525');
brackets(2).kl = rat('-9/10000000000');
brackets(2).kr = rat('-8/10000000000');
x = iv(XLEFT, XRIGHT, C);
y = iv(YLOW, YHIGH, C);
z = iv(Z, C);
fprintf('Singular-kite certificates for Theorem 7.2 and the local mass surface\n');
fprintf('z_* = %s\n', ratText(Z));
fprintf('rectangle x %s, y %s\n', ivText(x, C), ivText(y, C));
fprintf('subdivision = %d x %d exact rational boxes\n\n', ...
        C.subdivisions, C.subdivisions);
% Certify the complete rectangle by an exact rational subdivision.  The
% union of the subboxes is exactly [XLEFT,XRIGHT] x [YLOW,YHIGH].
bounds = certifyRectangle(XLEFT, XRIGHT, YLOW, YHIGH, Z, ...
                          sqrt3, C.subdivisions, C);
% Membership in the semialgebraic branch Omega_{>=}^J.
% 3*(3/5)^2>1 proves 3/5>1/sqrt(3), while
% 2-SQRT3LOW < 67/250 proves 2-sqrt(3)<67/250.
cutoff = rat('67/250');
q = rat('3/5');
require(ratLT(rat('1'), ratMul(rat('3'), ratMul(q, q))));
require(ratLT(ratSub('2', SQRT3LOW), cutoff));
require(iinside(x, q, rat('1'), C));
require(iinside(y, cutoff, rat('1'), C));
require(iinside(z, rat('0'), rat('1'), C));
require(iinside(bounds.A, rat('0.9'), rat('0.91'), C));
require(iinside(bounds.B, rat('0.65'), rat('0.66'), C));
require(iinside(bounds.F, rat('1.159'), rat('1.161'), C));
require(iinside(bounds.H, rat('0.0075'), rat('0.0077'), C));
require(iinside(bounds.J, rat('1.82'), rat('1.83'), C));
require(ipos(isub(bounds.alpha, iv('1', C), C)));
fprintf('A = %s\n', ivText(bounds.A, C));
fprintf('B = %s\n', ivText(bounds.B, C));
fprintf('F = %s\n', ivText(bounds.F, C));
fprintf('H = %s\n', ivText(bounds.H, C));
fprintf('J = %s\n', ivText(bounds.J, C));
% Strict monotonicity and the regularity minor on the certified z=z_*
% slice.  The minor need only be nonzero at the common root for the
% analytic implicit-function argument.
require(iinside(bounds.Ex, rat('-2.455'), rat('-2.454'), C));
require(iinside(bounds.Ey, rat('-4.860'), rat('-4.858'), C));
require(iinside(bounds.Jac, rat('-1/25'), rat('-9/250'), C));
fprintf('E_x = %s\n', ivText(bounds.Ex, C));
fprintf('E_y = %s\n', ivText(bounds.Ey, C));
fprintf('det D_(x,y)(E,K3) = %s\n', ivText(bounds.Jac, C));
% Since E_x<0, the two corner signs control the entire horizontal faces:
% E(x,YLOW)>=E(XRIGHT,YLOW)>0 and
% E(x,YHIGH)<=E(XLEFT,YHIGH)<0 for every x in [XLEFT,XRIGHT].
bottom = kite(iv(XRIGHT, C), iv(YLOW, C), z, C);
top = kite(iv(XLEFT, C), iv(YHIGH, C), z, C);
bottomMin = bottom.E;
topMax = top.E;
require(ipos(bottomMin) && ineg(topMax));
require(iinside(bottomMin, rat('0.000000495'), rat('0.000000497'), C));
require(iinside(topMax, rat('-0.000000472'), rat('-0.000000469'), C));
fprintf('E(x_right,y_low) = %s\n', ivText(bottomMin, C));
fprintf('E(x_left,y_high) = %s\n', ivText(topMax, C));
% First IVT: for each x there is a zero y=h(x).  E_y<0 makes it unique;
% the analytic implicit function theorem makes h continuous (indeed
% analytic) locally, and uniqueness glues the local graphs across [a,b].
for k = 1:numel(brackets)
    xp = brackets(k).xp;
    yl = brackets(k).yl;
    yr = brackets(k).yr;
    kl = brackets(k).kl;
    kr = brackets(k).kr;
    require(ratLT(YLOW, yl) && ratLT(yl, yr) && ratLT(yr, YHIGH));
    leftEval = kite(iv(xp, C), iv(yl, C), z, C);
    rightEval = kite(iv(xp, C), iv(yr, C), z, C);
    el = leftEval.E;
    er = rightEval.E;
    require(ipos(el) && ineg(er));
    boxEval = kite(iv(xp, C), iv(yl, yr, C), z, C);
    Kbox = boxEval.K;
    require(iinside(Kbox, kl, kr, C));
    fprintf('endpoint x = %s, h(x) in %s\n', ...
            ratText(xp), ivText(iv(yl, yr, C), C));
    fprintf('  E at lower/upper y = %s  %s\n', ...
            ivText(el, C), ivText(er, C));
    fprintf('  K3 enclosure = %s\n', ivText(Kbox, C));
end
% Second IVT: K3(x,h(x),z_*) is continuous and has opposite signs at
% x=XLEFT and x=XRIGHT.  Hence an exact common zero p_* of E and K3 exists.
% Positivity of the singular mass family at p_*.  On the entire rectangle
% Delta2<0 and Delta5>0, so the signs of the two displayed numerators are
% exactly the signs required in the manuscript's parametrization.
require(ineg(bounds.d2) && ipos(bounds.d5));
require(iinside(bounds.d1, rat('-3'), rat('0'), C));
require(ipos(isub(inegative(iv('1/4', C)), bounds.U, C)));
require(ipos(isub(bounds.W, iv('1/20', C), C)));
require(iinside(bounds.T, rat('0'), rat('4'), C));
fprintf(['coarse positivity bounds: -3<Delta1<0, U<-1/4, ' ...
         'W>1/20, 0<T<4\n']);
fprintf('Delta1 = %s\n', ivText(bounds.d1, C));
fprintf('U = %s\n', ivText(bounds.U, C));
fprintf('W = %s\n', ivText(bounds.W, C));
fprintf('T = (alpha-1)Delta6 = %s\n', ivText(bounds.T, C));
% For every 0<t<=1/100,
% U-Delta1*t < -1/4+3/100 = -11/50 < 0,
% W-T*t > 1/20-4/100 = 1/100 > 0,
% and m3=t>0.  Thus the whole interval of masses is positive.
whole = kite(x, y, z, C);
t = iv('1/100', C);
[m1, m5] = singularMasses(whole, t, C);
require(iinside(m1, rat('0.0172'), rat('0.0173'), C));
require(iinside(m5, rat('0.3018'), rat('0.3020'), C));
fprintf('m1 at t=1/100 = %s\n', ivText(m1, C));
fprintf('m5 at t=1/100 = %s\n', ivText(m5, C));
% These additional bounds use the entire fixed-z_* rectangle, including
% its exact common zero, rather than a rounded approximation of the root.
extensions = certifyExtensions(x, y, z, C);
fprintf('\nFull positive-mass interval and its endpoints:\n');
fprintf('tau = W/((alpha-1)Delta6) = %s\n', ivText(extensions.tau, C));
fprintf('rho = U/Delta1 = %s\n', ivText(extensions.rho, C));
fprintf('at t=0: m3=0, m1 in %s, m5 in %s\n', ...
        ivText(extensions.m1AtZero, C), ivText(extensions.m5AtZero, C));
fprintf('at t=tau: m1=0 exactly by definition, m3=tau>0, m5 in %s\n', ...
        ivText(extensions.m5AtTau, C));
fprintf('\nMass-surface derivative certificate at t_*=1/100:\n');
fprintf('whole-rectangle det D_(x,y)(E,K3) = %s\n', ...
        ivText(extensions.Jac, C));
fprintf('x''(z) enclosure = %s\n', ivText(extensions.vx, C));
fprintf('y''(z) enclosure = %s\n', ivText(extensions.vy, C));
fprintf('partial_z S1 along sigma(z), with t fixed = %s\n', ...
        ivText(extensions.dS1dz, C));
fprintf('\nPASS: two one-dimensional IVTs prove an exact common zero E=K3=0.\n');
fprintf('PASS: the common zero lies in the manuscript domain Omega_adm^J.\n');
fprintf('PASS: I(p_*)=(0,tau(p_*)), with 0.0155<tau(p_*)<0.0157.\n');
fprintf('PASS: every 0<t<=1/100 gives a strictly positive mass vector.\n');
fprintf('PASS: the interval endpoints have m3=0 and m1=0; all other masses are positive.\n');
fprintf(['PASS: the nonzero Jacobian minor at the root gives the local ' ...
         'analytic singular-compatible curve.\n']);
fprintf('PASS: 179<partial_z S1(z_*,1/100)<197, with m3=t.\n');
fprintf(['PASS: the (m1,m3) projection is locally invertible; the singular ' ...
         'mass family is a local embedded real-analytic surface.\n']);
fprintf(['NOTE: the strict bounds extend along the nearby curve by ' ...
         'continuity after shrinking the analytic neighborhood.\n']);
fprintf('NOTE: no global uniqueness in the full rectangle is claimed.\n');
fprintf('NOTE: the mass-surface conclusion is local; no global injectivity is claimed.\n');
end

function out = certifyExtensions(x, y, z, C)
% Certify the maximal positive interval and the local mass projection.
% x,y cover the entire rectangle and z is the exact singleton z_*.
% The previously certified common zero lies in this rectangle.
zero = iv('0', C);
one = iv('1', C);
tstar = iv('1/100', C);
r = kite(x, y, z, C);
T = imul(isub(r.alpha, one, C), r.d{6}, C);
require(ineg(r.d{1}) && ineg(r.b{1}));
require(ipos(r.b{3}) && ipos(T));
require(ineg(r.d{2}) && ipos(r.d{5}) && ipos(r.alpha));
tau = idiv(r.b{3}, T, C);
rho = idiv(r.b{1}, r.d{1}, C);
require(iinside(tau, rat('0.0155'), rat('0.0157'), C));
require(iinside(rho, rat('0.095'), rat('0.097'), C));
require(ipos(isub(rho, tau, C)));
require(ipos(isub(tau, tstar, C)));
% For a compatible singular shape the signs above give precisely
% t>0, t<rho, t<tau.  Since tau<rho, I(p)=(0,tau(p)).
[m1AtZero, m5AtZero] = singularMasses(r, zero, C);
[~, m5AtTau] = singularMasses(r, tau, C);
require(ipos(m1AtZero) && ipos(m5AtZero));
require(ipos(tau) && ipos(m5AtTau));
% m1(p,tau(p))=0 is an exact algebraic identity.  Do not test this
% cancellation by requiring a dependency-inflated interval to equal zero.
verifyMassEndpointIdentity();
% Reuse the existing two-component dual arithmetic.  The first pass
% seeds x,y; the second seeds z alone.  All three partial derivatives are
% therefore available without changing the interval/dual primitives.
rxy = kite(dual(x, one, zero, C), dual(y, zero, one, C), z, C);
rz = kite(x, y, dual(z, one, zero, C), C);
Ex = rxy.E.d{1}; Ey = rxy.E.d{2};
Kx = rxy.K.d{1}; Ky = rxy.K.d{2};
Ez = rz.E.d{1}; Kz = rz.K.d{1};
Jac = isub(imul(Ex, Ky, C), imul(Ey, Kx, C), C);
require(iinside(Jac, rat('-1/25'), rat('-9/250'), C));
% D_(x,y)(E,K3)*(x',y')' = -(E_z,(K3)_z)'.
% These rational expressions equal the curve tangent at compatible roots.
vx = idiv(isub(imul(Ey, Kz, C), imul(Ez, Ky, C), C), Jac, C);
vy = idiv(isub(imul(Kx, Ez, C), imul(Ex, Kz, C), C), Jac, C);
[m1xy, ~] = singularMasses(rxy, tstar, C);
[m1z, ~] = singularMasses(rz, tstar, C);
dS1dz = iadd(iadd(imul(m1xy.d{1}, vx, C), ...
                  imul(m1xy.d{2}, vy, C)), m1z.d{1});
require(iinside(dS1dz, rat('179'), rat('197'), C));
% S2(z,t)=m3=t, so det D_(z,t)(S1,S2)=partial_z S1.
% The inverse function theorem turns this nonzero minor into a local
% graph S3=g(S1,S2), hence a local embedded real-analytic mass surface.
out = struct();
out.tau = tau;
out.rho = rho;
out.m1AtZero = m1AtZero;
out.m5AtZero = m5AtZero;
out.m5AtTau = m5AtTau;
out.Jac = Jac;
out.vx = vx;
out.vy = vy;
out.dS1dz = dS1dz;
end

function [m1, m5] = singularMasses(r, t, C)
% Evaluate (7.2) over ordinary or dual intervals; m3=t identically.
T = gmul(gsub(r.alpha, '1', C), r.d{6}, C);
m1 = gdiv(gsub(r.b{3}, gmul(T, t, C), C), ...
          gmul(r.alpha, r.d{5}, C), C);
m5 = gdiv(gsub(r.b{1}, gmul(r.d{1}, t, C), C), ...
          gmul(r.alpha, r.d{2}, C), C);
end

function bounds = certifyRectangle(xlo, xhi, ylo, yhi, zrat, ...
                                   sqrt3, count, C)
% Cover the rectangle by exact rational subboxes and hull all enclosures.
xparts = ratPartition(xlo, xhi, count);
yparts = ratPartition(ylo, yhi, count);
first = true;
for ix = 1:count
    xb = iv(xparts(ix).lo, xparts(ix).hi, C);
    for iy = 1:count
        yb = iv(yparts(iy).lo, yparts(iy).hi, C);
        zb = iv(zrat, C);
        r = kite(xb, yb, zb, C);
        dx = dual(xb, iv('1', C), iv('0', C), C);
        dy = dual(yb, iv('0', C), iv('1', C), C);
        rd = kite(dx, dy, zb, C);
        Ex = rd.E.d{1};
        Ey = rd.E.d{2};
        Kx = rd.K.d{1};
        Ky = rd.K.d{2};
        Jac = isub(imul(Ex, Ky, C), imul(Ey, Kx, C), C);
        F = gmul(gsub('1', gmul(yb, yb, C), C), ...
                 gadd('1', gmul(xb, zb, C), C), C);
        F = gsub(F, gmul(gmul('2', yb, C), ...
                           gsub(xb, zb, C), C), C);
        H = gsub(gsub(xb, zb, C), ...
                 gmul(gsub('2', sqrt3, C), ...
                      gadd('1', gmul(xb, zb, C), C), C), C);
        J = gmul(gsub('1', gmul(xb, xb, C), C), ...
                 gsub('1', gmul(zb, zb, C), C), C);
        J = gadd(J, gmul(gmul('4', xb, C), zb, C), C);
        J = gsub(J, gmul(gmul(gmul('2', yb, C), ...
                                    gsub(xb, zb, C), C), ...
                              gadd('1', gmul(xb, zb, C), C), C), C);
        T = imul(isub(r.alpha, iv('1', C), C), r.d{6}, C);
        values = {r.A, r.B, r.alpha, F, H, J, r.d{1}, r.d{2}, ...
                  r.d{5}, r.b{1}, r.b{3}, T, Ex, Ey, Jac};
        names = {'A', 'B', 'alpha', 'F', 'H', 'J', 'd1', 'd2', ...
                 'd5', 'U', 'W', 'T', 'Ex', 'Ey', 'Jac'};
        if first
            for k = 1:numel(names)
                bounds.(names{k}) = values{k};
            end
            first = false;
        else
            for k = 1:numel(names)
                bounds.(names{k}) = ihull(bounds.(names{k}), values{k});
            end
        end
    end
end
end

function parts = ratPartition(lo, hi, count)
% Exact rational partition; adjacent subintervals share exact endpoints.
step = ratDiv(ratSub(hi, lo), rat(sprintf('%d', count)));
for k = 1:count
    leftMultiplier = rat(sprintf('%d', k - 1));
    rightMultiplier = rat(sprintf('%d', k));
    parts(k).lo = ratAdd(lo, ratMul(leftMultiplier, step)); %#ok<AGROW>
    parts(k).hi = ratAdd(lo, ratMul(rightMultiplier, step)); %#ok<AGROW>
end
require(ratLE(parts(1).lo, lo) && ratLE(lo, parts(1).lo));
require(ratLE(parts(end).hi, hi) && ratLE(hi, parts(end).hi));
end

function out = ihull(a,b)
% Exact hull of two fixed-point interval enclosures.
out=ivRaw(min([a.lo,b.lo]),max([a.hi,b.hi]));
end

function out = kite(x, y, z, C)
% Evaluate the kite formula over intervals or two-component dual intervals.
A = gmul(gadd(x, y, C), gsub('1', gmul(x, y, C), C), C);
B = gmul(gsub(x, z, C), gadd('1', gmul(x, z, C), C), C);
onePlusXX = gadd('1', gmul(x, x, C), C);
onePlusYY = gadd('1', gmul(y, y, C), C);
onePlusZZ = gadd('1', gmul(z, z, C), C);
R = gmul(gmul(y, z, C), gpow(onePlusXX, 2, C), C);
alpha = gadd('1', gdiv(R, gmul(A, B, C), C), C);
beta = gdiv(gmul(gmul('2', y, C), ...
                     gsub('1', gmul(x, x, C), C), C), A, C);
gamma = gdiv(gmul(gmul('2', x, C), ...
                      gsub('1', gmul(y, y, C), C), C), A, C);
s12 = gpow(gdiv(A, gmul(x, onePlusYY, C), C), 3, C);
s23 = gpow(gdiv(A, gmul(y, onePlusXX, C), C), 3, C);
s24 = gpow(gdiv(A, gmul(gmul('4', x, C), y, C), C), 3, C);
den25 = gmul(gmul(gmul(x, y, C), onePlusXX, C), onePlusZZ, C);
s25 = gpow(gdiv(gmul(A, B, C), den25, C), 3, C);
s35 = gpow(gdiv(gmul(A, B, C), R, C), 3, C);
s15 = gpow(alpha, -3, C);
d1 = gsub('1', s23, C);
d2 = gsub(s15, s25, C);
d3 = gsub(s12, '1', C);
d4 = gsub(s25, s35, C);
d5 = gsub(s12, s15, C);
d6 = gsub(s23, s35, C);
U = gmul(gsub(s24, s12, C), gamma, C);
V = gmul(gsub(s23, s24, C), beta, C);
W = gmul(gsub(s25, s24, C), ...
         gsub(gadd(beta, gmul('2', alpha, C), C), '2', C), C);
E = gsub(gmul(gmul(d1, d4, C), d5, C), ...
         gmul(gmul(d2, d3, C), d6, C), C);
Kterm1 = gmul(gmul(d2, d3, C), W, C);
Kterm2 = gmul(gmul(gmul(gsub('1', alpha, C), d4, C), d5, C), U, C);
Kterm3 = gmul(gmul(gmul(alpha, d2, C), d5, C), V, C);
K = gsub(gadd(Kterm1, Kterm2, C), Kterm3, C);
out.E = E;
out.K = K;
out.alpha = alpha;
out.d = {d1, d2, d3, d4, d5, d6};
out.b = {U, V, W};
out.A = A;
out.B = B;
end
% ---------- Generic interval/dual arithmetic ----------

function out = gadd(a, b, C)
if isDual(a) || isDual(b)
    a = asDual(a, C); b = asDual(b, C);
    out = dual(iadd(a.v, b.v), ...
               iadd(a.d{1}, b.d{1}), iadd(a.d{2}, b.d{2}), C);
else
    out = iadd(iv(a, C), iv(b, C));
end
end

function out = gnegative(a, C)
if isDual(a)
    out = dual(inegative(a.v), inegative(a.d{1}), ...
               inegative(a.d{2}), C);
else
    out = inegative(iv(a, C));
end
end

function out = gsub(a, b, C)
out = gadd(a, gnegative(b, C), C);
end

function out = gmul(a, b, C)
if isDual(a) || isDual(b)
    a = asDual(a, C); b = asDual(b, C);
    d1 = iadd(imul(a.d{1}, b.v, C), imul(b.d{1}, a.v, C));
    d2 = iadd(imul(a.d{2}, b.v, C), imul(b.d{2}, a.v, C));
    out = dual(imul(a.v, b.v, C), d1, d2, C);
else
    out = imul(iv(a, C), iv(b, C), C);
end
end

function out = gdiv(a, b, C)
if isDual(a) || isDual(b)
    a = asDual(a, C); b = asDual(b, C);
    den = ipow(b.v, 2, C);
    d1 = idiv(isub(imul(a.d{1}, b.v, C), ...
                      imul(b.d{1}, a.v, C), C), den, C);
    d2 = idiv(isub(imul(a.d{2}, b.v, C), ...
                      imul(b.d{2}, a.v, C), C), den, C);
    out = dual(idiv(a.v, b.v, C), d1, d2, C);
else
    out = idiv(iv(a, C), iv(b, C), C);
end
end

function out = gpow(a, n, C)
require(isscalar(n) && n == fix(n));
if isDual(a)
    if n == 0
        out = dual(iv('1', C), iv('0', C), iv('0', C), C);
        return
    end
    factor = imul(iv(sprintf('%d', n), C), ipow(a.v, n - 1, C), C);
    out = dual(ipow(a.v, n, C), imul(factor, a.d{1}, C), ...
               imul(factor, a.d{2}, C), C);
else
    out = ipow(iv(a, C), n, C);
end
end

function tf = isDual(x)
tf = isstruct(x) && isfield(x, 'kind') && strcmp(x.kind, 'dual');
end

function out = asDual(x, C)
if isDual(x)
    out = x;
else
    out = dual(iv(x, C), iv('0', C), iv('0', C), C);
end
end

function out = dual(value, dx, dy, C)
out = struct();
out.kind = 'dual';
out.v = iv(value, C);
out.d = {iv(dx, C), iv(dy, C)};
end
% ---------- Symbolic fixed-point outward-rounded intervals ----------

function out = iv(lower,upper,C)
% Convert exact rational bounds to a 60-place outward fixed-point interval.
if nargin==2
    C=upper;
    upper=lower;
end
if isstruct(lower) && isfield(lower,'kind') && strcmp(lower.kind,'interval')
    require(nargin==2,'Cannot give a second endpoint when copying an interval.');
    out=lower;
    return
end
loRat=rat(lower); hiRat=rat(upper);
require(loRat<=hiRat,'Interval endpoints are reversed.');
lo=floor(loRat*C.scale);
hi=ceil(hiRat*C.scale);
out=ivRaw(lo,hi);
end

function out = ivRaw(lo,hi)
% lo and hi are exact symbolic integers, representing endpoints / C.scale.
out=struct('kind','interval','lo',lo,'hi',hi);
end

function out = iadd(a,b)
out=ivRaw(a.lo+b.lo,a.hi+b.hi);
end

function out = inegative(a)
out=ivRaw(-a.hi,-a.lo);
end

function out = isub(a, b, ~)
out = iadd(a, inegative(b));
end

function out = imul(a,b,C)
p=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi];
out=ivRaw(floorDiv(min(p),C.scale),ceilDiv(max(p),C.scale));
end

function out = ireciprocal(a,C)
require((a.lo>0) | (a.hi<0),'Interval denominator contains zero.');
out=ivRaw(floorDiv(C.scale2,a.hi),ceilDiv(C.scale2,a.lo));
end

function out = idiv(a, b, C)
out = imul(a, ireciprocal(b, C), C);
end

function out = ipow(a, n, C)
require(isscalar(n) && n == fix(n));
if n < 0
    out = ipow(ireciprocal(a, C), -n, C);
    return
end
out = iv('1', C);
for k = 1:n
    out = imul(out, a, C);
end
end

function tf = iinside(a,lower,upper,C)
lower=rat(lower); upper=rat(upper);
tf=proved((a.lo>lower*C.scale) & (a.hi<upper*C.scale));
end

function tf = ipos(a)
tf=proved(a.lo>0);
end

function tf = ineg(a)
tf=proved(a.hi<0);
end

function out = rat(value)
% Local exact parser; this is not MATLAB's floating-point rat approximation.
if isa(value,'sym')
    require(isscalar(value),'Rational input must be scalar.');
    out=value;
    return
end
if isnumeric(value) || islogical(value)
    require(isscalar(value) && isreal(value),'Numeric rational input must be scalar.');
    if isinteger(value)
        out=sym(char(string(value)));
    else
        require(isfinite(value) && value==fix(value) && abs(value)<=flintmax, ...
            'Numeric rational inputs must be safe integers; use a string otherwise.');
        out=sym(value);
    end
    return
end
if isstring(value)
    require(isscalar(value),'Expected a scalar string.');
    value=char(value);
end
require(ischar(value) && isrow(value),'Expected a rational/decimal string.');
s=strtrim(value);
require(~isempty(s),'Empty rational string.');
slash=strfind(s,'/');
if ~isempty(slash)
    require(numel(slash)==1,'Invalid rational string.');
    ns=strtrim(s(1:slash-1)); ds=strtrim(s(slash+1:end));
    require(~isempty(regexp(ns,'^[+-]?[0-9]+$','once')) && ...
        ~isempty(regexp(ds,'^[+-]?[0-9]+$','once')),'Invalid fraction.');
    denominator=sym(ds); require(denominator~=0,'Zero rational denominator.');
    out=sym(ns)/denominator;
else
    require(~isempty(regexp(s,'^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)$','once')), ...
        'Invalid decimal rational.');
    dot=strfind(s,'.');
    if isempty(dot)
        out=sym(s);
    else
        places=numel(s)-dot; s(dot)=[];
        out=sym(s)/sym(10)^places;
    end
end
end

function out = ratMul(a,b)
out=rat(a)*rat(b);
end

function out = ratAdd(a,b)
out=rat(a)+rat(b);
end

function out = ratSub(a,b)
out=rat(a)-rat(b);
end

function out = ratDiv(a,b)
a=rat(a); b=rat(b);
require(b~=0,'Zero rational divisor.');
out=a/b;
end

function tf = ratLT(a,b)
tf=proved(rat(a)<rat(b));
end

function tf = ratLE(a,b)
tf=proved(rat(a)<=rat(b));
end

function q = floorDiv(a,b)
% Mathematical floor of an exact integer ratio, including negative inputs.
require(b~=0,'Zero integer denominator.');
q=floor(a/b);
end

function q = ceilDiv(a,b)
require(b~=0,'Zero integer denominator.');
q=ceil(a/b);
end

function s = ratText(value)
[n,d]=numden(rat(value));
if proved(d==1), s=char(n);
else, s=[char(n) '/' char(d)]; end
end

function s = ivText(a, C)
s = ['[' scaledText(a.lo, C) ', ' scaledText(a.hi, C) ...
     '] (display rounded)'];
end

function s = scaledText(x,C)
% Round only the printed value to 18 digits after the decimal point.
shown=18; drop=C.prec-shown;
require(drop>=0,'Insufficient internal decimal precision.');
factor=sym(10)^drop;
ax=abs(x); q=floor(ax/factor); r=ax-q*factor;
if proved(2*r>=factor), q=q+1; end
digits=char(q);
if numel(digits)<=shown
    digits=[repmat('0',1,shown+1-numel(digits)) digits];
end
intPart=digits(1:end-shown); fracPart=digits(end-shown+1:end);
while ~isempty(fracPart) && fracPart(end)=='0', fracPart(end)=[]; end
if isempty(fracPart), s=intPart;
else, s=[intPart '.' fracPart]; end
if proved(x<0) && proved(q~=0), s=['-' s]; end
end

function verifyMassEndpointIdentity()
% Symbolically prove the endpoint cancellation used in the interval argument.
% Positivity checks in certifyExtensions certify the actual denominators.
syms alpha0 Delta5 W T t real
m1=(W-T*t)/(alpha0*Delta5);
equalZero(subs(m1,t,W/T),'The exact identity m1(p,tau(p))=0');
end

function require(condition,label)
% Accept only conditions that Symbolic Math Toolbox can prove true.
if nargin<2, label='An exact certificate check failed.'; end
if isa(condition,'sym')
    truth=false(size(condition));
    for k=1:numel(condition)
        truth(k)=isAlways(condition(k),'Unknown','false');
    end
else
    truth=logical(condition);
end
if ~all(truth(:))
    error('certify_singular_kite_ivt:CheckFailed','%s',label);
end
end

function tf = proved(condition)
if isa(condition,'sym')
    tf=all(isAlways(condition(:),'Unknown','false'));
else
    tf=all(logical(condition(:)));
end
end

function equalZero(expression,label)
% Prove each rational expression is zero via its numerator.
for k=1:numel(expression)
    reduced=simplifyFraction(expression(k));
    [numerator,~]=numden(reduced);
    numerator=simplify(expand(numerator));
    require(numerator==0,label);
end
end
