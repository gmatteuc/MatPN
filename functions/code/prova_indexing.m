% %% folders
% 
% listSessions = {};
% expected_num = 13;
% num = 0;
% 
% rootFolder = fullfile('/zocconasphys1','chronic_inv_rec','Tanks');
% list = dir(rootFolder);
% for i=1:numel(list)
%    if( strfind(list(i).name, 'Giulio_Acute_Recording') ) 
%             
%       parts = strsplit(list(i).name, '_');
%       
%       if ( numel(parts) == 6 & strcmp('2016', parts{end}) )
%           num = num + 1;         
%           listSessions{num} = fullfile(rootFolder, list(i).name); 
%          
%       end
%       
%    end
% end
% 
% warning('expected = %d --- found = %d', expected_num, num)
% if(num ~= expected_num)   
%    error('Wrong number of sessions !') 
% else
%     
% end