
% RESET THE DIAGNOSTIC VARIABLES
%             if obj.objectID == 1
%                 obj.LP1 = NaN(2,1);
%                 obj.LP2 = NaN(2,1);
%                 obj.LP3 = NaN(2,1);
%                 obj.testVariable = NaN(2,1);
%             end

%             if obj.TIME.currentTime == 9.5
%                display('problem time');
%                obj.orcaLines
%             end

% ////////// WRITE DATA TO CSV //////////////////////////
%             if obj.objectID == 1
%                 params = [obj.TIME.currentTime,obj.objectID,...
%                           obj.VIRTUAL.globalPosition(1),...
%                           obj.VIRTUAL.globalPosition(2),...
%                           obj.VIRTUAL.globalVelocity(1),...
%                           obj.VIRTUAL.globalVelocity(2),...
%                           8888,...
%                           1000,obj.LP1(1),obj.LP1(2),...
%                           2000,obj.LP2(1),obj.LP2(2),...
%                           3000,obj.LP3(1),obj.LP3(2)];
% %                           8888,desiredVelocity(1),desiredVelocity(2),...
% %                           8888,obj.testVariable(1)];
%                 % WRITE TO CSV FILE
%                 obj.sendVectorToCsv(obj.TIME,params)
%             end