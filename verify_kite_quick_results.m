function verify_kite_quick_results()
%VERIFY_KITE_QUICK_RESULTS Rigorous certificate for Theorem 9.1 and 9.2.
%
% Basis: the exterior-fifth-body parametrization and linear system
% (3.6)-(3.8).
%
% Requires Symbolic Math Toolbox. Run from MATLAB with:
%     verify_kite_quick_results
% or run all bundled certificates with:
%     run_all
%
% Rational arithmetic, not floating-point tolerances, certifies all
% inequalities. Decimal values are printed only for readability. This is a
% local certificate; it does NOT certify global connectedness or
% nonexistence of singular positive solutions. The accompanying geometric
% proofs remain necessary.

if isempty(ver('symbolic'))
    error('verify_kite_quick_results:MissingSymbolicToolbox', ...
          'This verification script requires Symbolic Math Toolbox.');
end

syms x y z real
vars = [x, y, z];
A = (x + y)*(1 - x*y);
B = (x - z)*(1 + x*z);
R = y*z*(1 + x*x)^2;
C = A*B + R;
F = (1 - y*y)*(1 + x*z) - 2*y*(x - z);
G = x*(1 - z*z) - 2*z;
J = (1 - x*x)*(1 - z*z) + 4*x*z - 2*y*B;

% The new necessary inequality J>0 strictly strengthens the declared
% domain.
root3 = sqrt(sym(3));
alpha = C/(A*B);
beta = 2*y*(1 - x*x)/A;
gamma = 2*x*(1 - y*y)/A;
r12 = x*(1 + y*y)/A;
r23 = y*(1 + x*x)/A;
r24 = 4*x*y/A;
r25 = x*y*(1 + x*x)*(1 + z*z)/(A*B);
r35 = R/(A*B);
s12 = r12^(-3);
s23 = r23^(-3);
s24 = r24^(-3);
s25 = r25^(-3);
s35 = r35^(-3);
s15 = alpha^(-3);
d1 = 1 - s23;
d2 = s15 - s25;
d3 = s12 - 1;
d4 = s25 - s35;
d5 = s12 - s15;
d6 = s23 - s35;
P = [0,        d1,             alpha*d2;
     d3,       0,              (1 - alpha)*d4;
     alpha*d5, (alpha - 1)*d6, 0];
b = [(s24 - s12)*gamma;
     (s23 - s24)*beta;
     (s25 - s24)*(beta + 2*alpha - 2)];
equalZero(d5 - d3 - (1 - s15), 'Rank-cover identity');

p0 = [sym(1)/2, sym(1)/2, sym(1)/5];
require([2 - root3 < p0(1); p0(1) < 1/root3], 'x lower branch');
require([2 - root3 < p0(2); p0(2) < 1], 'y bounds');
require([0 < p0(3); p0(3) < 1], 'z bounds');
F0 = subs(F, vars, p0);
G0 = subs(G, vars, p0);
J0 = subs(J, vars, p0);
require([F0 > 0; G0 > 0; J0 > 0], 'Admissibility including J>0');

P0 = subs(P, vars, p0);
b0 = subs(b, vars, p0);
% MATLAB factor(q) returns a vector of prime factors when q is a rational
% constant. Keep the determinant as one exact rational scalar instead.
detP = simplifyFraction(det(P0));
require(isscalar(detP), 'det P must be scalar');
require(detP ~= 0, 'Nonsingularity');
m0 = simplify(P0\b0);
equalZero(P0*m0 - b0, 'Exact linear-system residual');
bounds = [sym(98)/100,  sym(1);
          sym(166)/100, sym(168)/100;
          sym(26)/100,  sym(27)/100];
for k = 1:3
    require([bounds(k,1) < m0(k); m0(k) < bounds(k,2)], ...
            'Exact positive-mass bounds');
end
fprintf('PASS: p0=(1/2,1/2,1/5) is strictly admissible and has positive masses.\n');
fprintf('F(p0), G(p0), J(p0) =\n');
disp([F0, G0, J0]);

fprintf('det P(p0) =\n');
disp(detP);
fprintf('Mass vector, exact =\n');
disp(m0.');
fprintf('Mass vector, decimals =\n');
disp(vpa(m0.', 15));

% Differentiate P(p)m(p)=b(p), avoiding a huge symbolic matrix inverse.
Dmass = sym(zeros(3,3));
for k = 1:3
    variable = vars(k);
    rhs = subs(diff(b, variable), vars, p0) ...
          - subs(diff(P, variable), vars, p0)*m0;
    Dmass(:,k) = simplify(P0\rhs);
end
detDmass = simplifyFraction(det(Dmass));
require(isscalar(detDmass), 'Mass Jacobian determinant must be scalar');
require([sym(2899) < detDmass; detDmass < sym(2900)], ...
        'Exact Jacobian determinant certificate');
fprintf('PASS: 2899 < det D(m1,m3,m5)(p0) < 2900 (exact rational comparison).\n');
fprintf('Mass Jacobian, decimal display =\n');
disp(vpa(Dmass, 12));
fprintf('Mass Jacobian determinant, exact =\n');
disp(detDmass);
fprintf('Mass Jacobian determinant, decimal =\n');
disp(vpa(detDmass, 18));

fprintf('All exact algebraic and local checks completed.\n');
end


function require(condition, label)
% Accept only conditions that Symbolic Math Toolbox can prove true.
if isa(condition, 'sym')
    truth = false(size(condition));
    for k = 1:numel(condition)
        truth(k) = isAlways(condition(k), 'Unknown', 'false');
    end
else
    truth = logical(condition);
end
if ~all(truth(:))
    error('verify_kite_quick_results:CheckFailed', '%s', label);
end
end


function equalZero(expression, label)
% Prove each rational expression is identically zero via its numerator.
for k = 1:numel(expression)
    reduced = simplifyFraction(expression(k));
    [numerator, ~] = numden(reduced);
    numerator = simplify(expand(numerator));
    require(numerator == 0, label);
end
end
