%% Laminate Class
% laminas = array of objects of Lamina Class
% status = specifies if lamina is failed or not
% z = distance from the midsurface

classdef Laminate
    properties
        laminas
        status
        z
    end
    methods
        function obj = Laminate(geo, elas, stren, hygro)
            if(size(geo,1)==size(elas,1)&&size(elas,1)==size(stren,1)&&size(stren,1)==size(hygro,1))
                obj.laminas = [];
                for ii = 1:size(geo, 1)
                    obj.laminas = [obj.laminas Lamina(geo(ii,:), elas(ii,:), stren(ii,:), hygro(ii,:))];
                end
                obj.status = ones(1, size(geo, 1));
                mid = sum(geo(:,1))/2;
                obj.z = zeros(1, length(obj.laminas) + 1);
                obj.z(1) = -mid;
                for ii = 2:length(obj.z)
                    mid = mid - obj.laminas(ii - 1).geometric(1);
                    obj.z(ii) = -mid;
                end
            else
                str_1 = "number of rows in input varibles are not consistant";
                str_2 = "check if geometric, elastic, strength and hygrothermal data have same number of rows";
                error(str_1 + str_2)
            end
        end
        function res = is_failed(obj)
            % returns true if all the lamina in the laminate have failed
            res = (sum(obj.status(obj.status ~= 0)) == 0);
        end
        function out = ABD(obj)
            % return ABD Matrix of the laminate
            A = zeros(3);
            B = zeros(3);
            D = zeros(3);
            for ii = 1:length(obj.laminas)
                A = A + obj.laminas(ii).Qbar*(obj.z(ii + 1) - obj.z(ii));
                B = B + 0.5*obj.laminas(ii).Qbar*(obj.z(ii + 1)^2 - obj.z(ii)^2);
                D = D + (1/3)*obj.laminas(ii).Qbar*(obj.z(ii + 1)^3 - obj.z(ii)^3);
            end
            out = [A, B; B, D];
        end
    end
end


%% First and Last Ply Failure Loads

% Coding Assignment
% ME607: Introduction to Composite Materials

% Name = Mayank Pathania
% Roll No. = 204103314
% Specialization = Machine Design
% Indian Institute of Technology, Guwahati
