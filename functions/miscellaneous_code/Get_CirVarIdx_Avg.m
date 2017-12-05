function [ Avg_R, Avg_T ] = Get_CirVarIdx_Avg( Vars, FR, Type )
%% GET CIRCULAR VARIANCE INFO ABOUT F.R. MATRIX FOR TYPE ('ORI','DIR')
%  GET R AND THETA FOR F.R. MATRIX FOR STIMULUS TIMES

%% --------------- CHECK IF FUNCTION IS CALLED FOR ORIENTATION OR DIRECTION
switch Type                                                                 % Check Type
    case 'Ori'                                                              % For Orientation
        Img = 2i;   Angle = 2;  FR = Get_FR_Dir2Ori( FR );                  % Set Variables and Transform FR 
    case 'Dir'                                                              % For Direction
       Img = 1i;    Angle = 1;                                              % Set Variables
end

%% -------------------------------------------------- SET UTILITY VARIABLES       
Radians = Vars.Radians;                                                     % Set Direction Vector in Radians
Dir_Tot = size( FR, 1 );                                                    % Get F.R. Matrix Directions (Rows)
Avg_Num = 0; Avg_Den = 0;                                                   % Reset CirVar Variables
if size(FR,2) > 1                                                           % Check if FR is a Matrix or a Vector
        FR_Avg_Dir = Get_FR_Avg_Dir( FR );                                  % Get Avg F.R. for each Direction
else    FR_Avg_Dir = FR;                                                    % Set FR with given Avg F.R. for each Direction
end
%% ----------------------------- GET AVERAGE CIRCULAR VARIANCE FOR STIMULUS
for Dir_Ind = 1 : Dir_Tot                                                   % For all Directions
    Avg_Num = Avg_Num + (FR_Avg_Dir(Dir_Ind) ...
        * exp( Img * Radians( Dir_Ind )));                                  % Get Avg CV Num.for specified direction
    Avg_Den = Avg_Den +  FR_Avg_Dir(Dir_Ind);                               % Get Avg CV Den. for specified direction
end
Avg_Temp = Avg_Num / Avg_Den;                                               % Get temporary Avg CirVarIdx
Avg_R = abs( Avg_Temp );                                                    % Get Avg CV Radius
Avg_T = angle( Avg_Temp ) / Angle;                                          % Get Avg CV Theta Angle in Radians

end

