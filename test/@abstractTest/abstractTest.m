classdef (Abstract) abstractTest

    properties (Hidden)
        Data_ %<- make this it's own class?
        keys %<- make this it's own class?
    end

    % properties (Abstract)
    %     G
    %     c
    %     A
    %     b
    % end
    % 
    % properties (Abstract, Dependent)
    %     n
    %     nG
    %     nC
    % end

    % methods (Abstract)
    %     obj = minSum(obj1,obj2);
    %     obj = linMap(obj1,obj2);
    %     obj = intersect(obj1,obj2);
    % end

    methods
        % set/get
        function obj = set.Data_(obj,Data)
            obj.Data_ = abstractTest.checkData(Data);
        end


        % operations
        obj = plus(obj1,obj2)
    end

    methods
        function Data_ = checkData(Data)
            Data_ = Data; % <- TODO: add dim check...
        end
    end


end