% test script for NPA with inner-product 
% find the maximum violation for CHSH inequality
clear all

%% generate NPA relaxation
io_config = {[2,2],[2,2],[]};
level = 1;
extra_monomials = struct('status',{},'as',{},'ao',{},'bs',{},'bo',{},'cs',{},'co',{});
isComplex = false;
npa = NPA(io_config, level, extra_monomials, isComplex);

%% formulate objective function and constraints
% generate measurement operators for Alice and Bob
A = GeneratePartyOps('A', 2, 1);
B = GeneratePartyOps('B', 2, 1);

% take product for A and B (computing bipartite measurements)
for x = 1:2
    for y = 1:2
        AB(x,y) = ProductOp(A(x), B(y));
    end
end

% map operators to variables
a = npa.Operator2Variable(A);
b = npa.Operator2Variable(B);
ab = npa.Operator2Variable(AB);

% formulate CHSH
chsh = 4*(ab(1,1) + ab(1,2) + ab(2,1) - ab(2,2)) - 4 * (a(1) + b(1)) + 2;

% here, we use YALMIP which only expresses any optimisation as
% minimisation, hence we need to include a minus sign for maximisation
obj = -chsh;

% the only constraint is due to the structure of the moment matrix
constraints = [npa.npa_constraints];

%% optimisation proper
sol = optimize(constraints, obj);
S = value(chsh);
fprintf('The maximum violation is given by %.4f.\n', S)