function bh_gyro_boat_startup()
    
    % Add some folders to the MATLAB search path
    p = mfilename('fullpath');
    [folder,name,ext] = fileparts(p);

    addpath([folder],                        '-begin');
    addpath([folder,filesep,'THE_PICS'],     '-begin');
    addpath([folder,filesep,'THE_UTILITIES'],'-begin');

    % echo the first few elements of our search path
    sp             = path;
    TGT_SPLIT_CHAR = pathsep;
    C              = strsplit(sp, TGT_SPLIT_CHAR);
    
    fprintf('\n %s', repmat('*',1,50) );
    fprintf('\n Just added the following folders to the ');
    fprintf('\n HEAD of your search path: \n');
    fprintf('\n    ---> %s', C{1:3});
    fprintf('\n %s', repmat('*',1,50) );
    fprintf('\n ... we are finished HERE ---> %s\n',mfilename);
        
    % assert that we have a new enough version to run this demo
    % R2019b is the minimum release and corresponds to MATLAB 
    % version 9.7
    MIN_required_ML_version = '9.7';
    if(verLessThan('MATLAB',MIN_required_ML_version))
          % inform the user that he needs a NEWER release and then exit
          error('###_ERROR:  you need at least R2019b to run this demo');
    end
       
end
%_#########################################################################
% function LOC_assert_version(your_version_str)
%   %your_version_str   = version('-release');
%   is_too_old_version = verLessThan('matlab', 'R2016a');
%   if(is_too_old_version)
%       tmp_str = sprintf('Your MATLAB release is: <%s>\n But this DEMO needs at least <R2016a>',your_version_str);
%       errordlg(tmp_str,'RELEASE CHECK', 'modal');
%       % throw an error
%       error('###_ERROR:  you need at least R2016a to run this demo');
%   end
% end
%--------------------------------------------------------------------------