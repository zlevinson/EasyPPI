% ----------------------------------------------------------% 
% File name: PPI_input.m
% 
% Description:
% Object to handle PROLITH input parameters
% ----------------------------------------------------------%
classdef PPI_input
	properties
		value;
		id;
		func;
	end

	methods(Static)
		function [a, b] = check_inputs(a, b)
			if strcmp(class(a), 'PPI_input')
				a = a.value;
			end

			if strcmp(class(b), 'PPI_input')
				b = b.value;
			end
		end
	end

	methods
		function self = PPI_input(id, func)
			% Class constructor
			% id - int
			% PROLITH input ID
			%
			% func - function handle
			% Function used to get the current value of the parameter
			self.id = id;
			self.func = func;
		end

		% Getter function for the value property
		function value = get.value(self)
			value = self.func();
		end

		% Operator overloading for all operators that make sense.
		% We know that self.value will always be a scalar.
		function v = plus(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a+b;
		end

		function v = minus(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a-b;
		end

		function v = uminus(a)
			v = -1*a.value;
		end

		function v = uplus(a)
			v = a.value;
		end

		function v = times(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a*b;
		end

		function v = mtimes(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a*b;
		end

		function v = rdivide(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a / b;
		end

		function v = mrdivide(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a / b;
		end

		function v = power(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a^b;
		end

		function v = mpower(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a^b;
		end

		function v = lt(a, b)	
			[a, b] = PPI_input.check_inputs(a, b);
			v = a<b;
		end

		function v = gt(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a>b;
		end

		function v = le(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a<=b;
		end

		function v = ge(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a>=b;
		end

		function v = ne(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a~=b;
		end

		function v = eq(a, b)
			[a, b] = PPI_input.check_inputs(a, b);
			v = a==b;
		end
	end
end