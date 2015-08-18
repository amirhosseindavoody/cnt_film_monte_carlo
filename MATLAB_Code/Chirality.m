classdef Chirality
    %Chirality object
    properties
        n
        m
    end
    
    methods
        function obj = Chirality(r,c) 
         if nargin == 2
            obj(r,c) = Chirality;
         end
        end
    end
end