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
			if strcmp(class(a), 'PPI_input')
				v = a.value + b;
			else
				v = a + b.value;
			end
		end

		function v = minus(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value - b;
			else
				v = a - b.value;
			end
		end

		function v = uminus(a)
			v = -1*a.value;
		end

		function v = uplus(a)
			v = a.value;
		end

		function v = times(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value * b;
			else
				v = a * b.value;
			end
		end

		function v = mtimes(a, b)
			v = self.times(a, b);
		end

		function v = rdivide(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value / b;
			else
				v = a / b.value;
			end
		end

		function v = mrdivide(a, b)
			v = self.rdivide(a, b);
		end

		function v = power(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value ^ b;
			else
				v = a ^ b.value;
			end
		end

		function v = mpower(a, b)
			v = self.power(a, b);
		end

		function v = lt(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value < b;
			else
				v = a < b.value;
			end
		end

		function v = gt(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value > b;
			else
				v = a > b.value;
			end
		end

		function v = le(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value <= b;
			else
				v = a <= b.value;
			end
		end

		function v = ge(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value >= b;
			else
				v = a >= b.value;
			end
		end

		function v = ne(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value ~= b;
			else
				v = a ~= b.value;
			end
		end

		function v = eq(a, b)
			if strcmp(class(a), 'PPI_input')
				v = a.value == b;
			else
				v = a == b.value;
			end
		end
	end
end