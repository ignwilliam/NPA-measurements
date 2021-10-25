% test script for NPA with inner-product 
% find the maximum guessing probability given the CHSH violation
clear all
N = 15;

%% generate NPA relaxation
io_config = {[2,2],[2,2],[2]};
level = 2;
extra_monomials = struct('status',{},'as',{},'ao',{},'bs',{},'bo',{},'cs',{},'co',{});
isComplex = false;
npa = NPA(io_config, level, extra_monomials, isComplex);

%% formulate objective function and constraints
% generate measurement operators for Alice, Bob and Eve
A = GeneratePartyOps('A', 2, 1);
B = GeneratePartyOps('B', 2, 1);
C = GeneratePartyOps('C', 1, 1);

% take product for A and B (computing bipartite measurements)
for x = 1:2
    for y = 1:2
        AB(x,y) = ProductOp(A(x), B(y));
    end
end
AC = ProductOp(A(1),C);

% map operators to variables
a = npa.Operator2Variable(A);
b = npa.Operator2Variable(B);
c = npa.Operator2Variable(C);
ab = npa.Operator2Variable(AB);
ac = npa.Operator2Variable(AC);

% formulate CHSH
chsh = 4*(ab(1,1) + ab(1,2) + ab(2,1) - ab(2,2)) - 4 * (a(1) + b(1)) + 2;

% guessing probability
pg = 2*ac - a(1) - c + 1;

% here, we use YALMIP which only expresses any optimisation as
% minimisation, hence we need to include a minus sign for maximisation
obj = -pg;

S = linspace(2,sqrt(8),N);
GP = zeros(1,N);

for i = 1:N
    % subject to 
    constraints = [npa.npa_constraints, chsh == S(i)];

    %% optimisation proper
    sol = optimize(constraints, obj);
    GP(i) = value(pg);
end

S_analytical = linspace(2,sqrt(8),100);
GP_analytical = 0.5 * real(1 + sqrt(2 - S_analytical.^2/4));
plot(S_analytical, GP_analytical, S, GP, 'o');