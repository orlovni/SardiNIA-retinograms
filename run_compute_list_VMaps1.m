%from = '_extracted_skeletonized_CLAHE_1k_2k.log';
%to = '_computed_sigs_1k_2k';

from = '_extracted_skeletonized_AREDS_Category_1_Other_NEW';
to = '_computed_sigs_AREDS_Cat1Other_';

from = ['../' from];
to   = ['../' to];

curDir = pwd;
from = fullfile( curDir,from );
to   = fullfile( curDir,to );

compute_list_VMaps( from,to );