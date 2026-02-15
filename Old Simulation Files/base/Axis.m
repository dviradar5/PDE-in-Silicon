% %% Axis 
% % Enumeration class for the Cartesian axes x, y, z.
% 
% %%% Description
% % The only three instances |Axis.x|, |Axis.y|, and |Axis.z| of |Axis| are used
% % extensively in MaxwellFDFD to represent and specify the x-, y-, and z-axis of
% % the Cartesian coordinate system.
% %
% % |Axis| is a subclass of |Enumerated|, so it inherits all the methods of
% % |Enumerated|.  See <Enumerated.html |Enumerated|> for more details.
% 
% %%% Instances
% % * |Axis.x|: instance representing the x-axis
% % * |Axis.y|: instance representing the y-axis
% % * |Axis.z|: instance representing the z-axis
% 
% %%% Methods
% % Below, |axis| is an instance of |Axis|.
% %
% % * |[p, q, r] = cycle(axis)|: |[p, q, r]| is a cyclic permutation of |[Axis.x,
% % Axis.y, Axis.z]| satisfying |r == axis|.
% 
% %%% Example
% %   % Test conversion to integers and strings.
% %   fprintf('# of instances of Axis: %d\n', Axis.count);
% %   for w = Axis.elems
% %       fprintf('The integer value of Axis.%s is %d\n', char(w), int(w));
% %   end
% %
% %   % Test cyclic permutation.
% %   [p q r] = cycle(Axis.y);
% %   fprintf('The cyclic permutation of [x, y, z] beginning with y is [%s, %s, %s]\n', char(r), char(p), char(q));
% 
% %%% See Also
% % <Enumerated.html |Enumerated|>
% 
% classdef Axis < Enumerated
% 	enumeration
% 		x('x')
% 		y('y')
% 		z('z')
% 	end
% 
% 	methods (Static)
% 		function elems = elems(ind)
% 			elems = [Axis.x, Axis.y, Axis.z];
% 			if nargin > 0  % ind
% 				elems = elems(ind);
% 			end
% 		end
% 		
% 		function count = count()
% 			count = length(Axis.elems);
% 		end
% 	end
% 	
% 	methods
% 		function [p q r] = cycle(this)
% 			r = this;
% 			p = Axis.elems(mod(int(this), Axis.count) + 1);
% 			q = Axis.elems(mod(int(this)+1, Axis.count) + 1);
% 		end
% 	end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Axis < Enumerated
    enumeration
        x('x')
        y('y')
        z('z')
    end

    methods (Static)
        function elems = elems(ind)
            persistent elem_list;
            if isempty(elem_list)
                elem_list = [Axis.x, Axis.y, Axis.z];
            end
            if nargin > 0  % if an index is provided
                elems = elem_list(ind);
            else
                elems = elem_list;
            end
        end

        function count = count()
            count = numel(Axis.elems());
        end
    end

    methods
        function [p, q, r] = cycle(this)
            r = this;
            elems = Axis.elems(); % Get the list of enumeration values
            idx = find(elems == this); % Find the index of current enumeration
        
        % Ensure cyclic permutation
        p = elems(mod(idx, Axis.count()) + 1);
        q = elems(mod(idx + 1, Axis.count()) + 1);
%             p = elems(mod(int32(this) - 1, Axis.count()) + 1);
%             q = elems(mod(int32(this), Axis.count()) + 1);
        end
    end
end

% [p, q, r] = cycle(Axis.y);
% disp([char(r), char(p), char(q)]);

