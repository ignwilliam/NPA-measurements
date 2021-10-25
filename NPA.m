classdef NPA
	properties
		level
		io_config
		monomials
		isComplex

		moment_matrix
		ref_table
		index_matrix

		variables
		Gamma
		npa_constraints
	end

	methods
		%% constructor for the NPA with inner-product class
		function sdp = NPA(io_config, level, extra_monomials, isComplex)
			% Gram matrix of the input states, input/output configuration, and relaxation level
			sdp.io_config = io_config;
			sdp.level = level;
			sdp.isComplex = isComplex;

			% generate monomials
			sdp.monomials = GenerateOps(io_config, level);
			sdp.monomials = [sdp.monomials, extra_monomials];

			% generate the moment matrix
			[sdp.moment_matrix, sdp.ref_table, sdp.index_matrix] = GenerateMomentMatrix(sdp.monomials);

			% calculate sizes
			num_mono = length(sdp.monomials);
			num_var = length(sdp.ref_table);

			% generate SDP variables
			if isComplex
				sdp.variables = sdpvar(num_var, 1,'full', 'complex');
				sdp.Gamma = sdpvar(num_mono, num_mono, 'full', 'complex');
			else
				sdp.variables = sdpvar(num_var,1, 'full');
				sdp.Gamma = sdpvar(num_mono, num_mono);
			end

			% generate NPA constraints
			sdp.npa_constraints = MomentMatrixConstraints(sdp);
		end

		%% impose NPA and inner-product constraints
		function constraints = MomentMatrixConstraints(sdp)
			% initialise constraints
			constraints = [];

			num_monomials = length(sdp.monomials);


            for k = 1:num_monomials
                for l = 1:num_monomials
                    row = k;
                    col = l;

                    switch sdp.index_matrix(k,l)
                        case 0
                            constraints = [constraints, sdp.Gamma(row,col) == 0];
                        case 1
                            constraints = [constraints, sdp.Gamma(row,col) == 1];
                        otherwise
                            constraints = [constraints, ...
                                           sdp.Gamma(row,col) == sdp.variables(sdp.index_matrix(k,l),1)];
                    end
                end
            end
		

			if sdp.isComplex
				constraints = [constraints, sdp.Gamma == sdp.Gamma'];
				constraints = [constraints, sdp.Gamma + sdp.Gamma' >= 0];
			else
				constraints = [constraints, sdp.Gamma >= 0];
			end
		end

		function vars = Operator2Variable(sdp, ops)
			
			% calculate dimension of the matrix of operators, ops
			dims = size(ops);
			
			% flatten ops
			ops = reshape(ops,[1,numel(ops)]);

			% convert to variables
			for i = 1:length(ops)
				vars(i) = sdp.variables(IndexOps(ops(i),sdp.ref_table),1);
			end

			% reshaping the output
			vars = reshape(vars,dims);
		end



	end

end