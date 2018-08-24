% This is for No-Tiling (e.g. number of tiles is 1) only.
%
% run_chop_up_sigMatrix_noClasses( '_CLAHEregistr_NoClasses_wholeTorso',.75,'hehe' );
% run_chop_up_sigMatrix_noClasses( '_CLAHEregistr_NoClasses_wholeTorso' );
%
function run_chop_up_sigMatrix_noClasses( filename,varargin ),

TRportion = 0.5;
if nargin > 1,TRportion = varargin{1};end

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% I could save the TE (test) data as a single file (default option) using 
% the same path as for TR (train) data, OR save TE data as separate files.
% This would happen when there is a specified 4th parameter, the path to
% the TE data
if nargin > 2, 
 TEpath = varargin{2}; 
 clear TEpath
 separate_exp_files = logical(1);
 chop_up_AgeUniform_sigMatrix_noClasses( filename,TRportion,1,separate_exp_files );
else,
 chop_up_AgeUniform_sigMatrix_noClasses( filename,TRportion,1 );
end


%%chop_up_sigMatrix_noClasses( filename,-2700,1 );
%chop_up_AgeUniform_sigMatrix_noClasses( filename,TRportion,1 );

% Form the train-name:
trainName = [filename '_train'];

% Here I am excluding image_ids from the saved vars:

load(trainName);
save( trainName,'category_names','dataset_name','image_paths','signature_labels','signature_matrix' );

end % eofunc
