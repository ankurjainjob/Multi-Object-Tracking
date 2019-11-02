classdef SymbolicExpressionComparator <  matlab.unittest.constraints.Comparator
    
    methods(Access=protected)
        function bool = supportsContainer(~, value)
            bool = isa(value, 'sym');
        end
        
        function bool = containerSatisfiedBy(comparator, actual, expected)
            bool = comparator.isAlwaysEqual(actual, expected) || ...
                comparator.isSimplifiedEqual(actual, expected);
        end
    end
    
    methods(Access=private)
        function bool = isAlwaysEqual(~, actual, expected)
            result = isAlways(actual == expected, 'Unknown','false');
            bool = all(result(:));
        end
        
        function bool = isSimplifiedEqual(constraint, actual, expected)
            simplifiedResult = simplify(actual - expected, 'Steps', 1000);
            bool = constraint.isAlwaysEqual(simplifiedResult, 0);
        end
    end
end