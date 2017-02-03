function [groups, isolated] = find_islands_1(mpc)
%FIND_ISLANDS_1  Finds islands in a network
%   GROUPS = FIND_ISLANDS_1(MPC)
%   [GROUPS, ISOLATED] = FIND_ISLANDS_1(MPC)
%
%   Returns the islands in a network. The return value GROUPS
%   is a cell array of vectors of the bus indices for each island.
%   The second and optional return value ISOLATED is a vector of
%   indices of isolated buses that have no connecting branches.

%% define named indices into data matrices
define_constants;

%% find islands
nb  = size(mpc.bus, 1);     %% number of buses
nl  = size(mpc.branch, 1);  %% number of branches

e2i = sparse(mpc.bus(:, BUS_I), ones(nb, 1), 1:nb, max(mpc.bus(:, BUS_I)), 1);
C_on = sparse(1:nl, e2i(mpc.branch(:, F_BUS)), -mpc.branch(:, BR_STATUS), nl, nb) + ...
       sparse(1:nl, e2i(mpc.branch(:, T_BUS)),  mpc.branch(:, BR_STATUS), nl, nb);

if nnz(C_on)
    Y = C_on'*C_on + 1e-6*speye(nb); % Bus admittance matrix made of branches with 1 Ohm resistance
    IB = zeros(nb,1); % Vector with island number for each bus
    ni = 0; % Number of islands
    IBX = find(IB == 0); % Vector with bus indices with unknown island number
    while ~isempty(IBX)
        k = IBX(1); % Take the first unclasified bus
        Y(k,k) = Y(k,k) + 1; % Ground bus k with a resistor of 1 Ohm
        I = zeros(nb,1); I(k) = 1; % Inject current of 1 A into bus k
        V = Y\I; % Calculate bus voltages
        Y(k,k) = Y(k,k) - 1; % Take out the resistor at bus k
        ind = V > 0.001; % Find buses whose voltage is not zero (buses in an island would have a voltage of 1 V, the others would have 0 V)
        if length(ind(ind)) > 1
            ni = ni + 1;
            IB(ind) = ni; % Buses in ind are island ni
        else
            IB(ind) = -1; % Isolated bus
        end
        IBX = find(IB == 0);
    end
    groups = cell(ni,1);
    for i = 1:ni
        groups{i} = find(IB == i)';
    end
    isolated = find(IB == -1)';
else
    groups = [];
    isolated = 1:nb;
end
