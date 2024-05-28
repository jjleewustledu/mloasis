classdef Oasis3Bold
    %% line1
    %  line2
    %  
    %  Created 05-Jan-2024 15:51:21 by jjlee in repository /Users/jjlee/MATLAB-Drive/mloasis/src/+mloasis.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2024 John J. Lee.
    
    properties
        sourcedata
    end

    methods
        function this = Oasis3Bold()
            this.sourcedata = "login3.chpc.wustl.edu:/scratch/jjlee/Singularity/OASIS3/bids2/sourcedata";
        end

        function mg = mglob_bold(~)
            %% N = 2138

            pwd0 = pushd("/data/oasis/OASIS3/MR/bold");
            mg = mglob("OAS3*_MR_d*");
            mg = strip(mg, filesep);
            popd(pwd0);    
        end
        function mgd = mglob_dwi(this)
            mgb = this.mglob_bold();
            pwd0 = pushd("/data/oasis/OASIS3/MR/dwi");
            mgd_ = mglob("OAS3*_MR_d*");
            mgd_ = strip(mgd_, filesep);
            j = 1;
            for i = 1:length(mgd_)
                re = regexp(mgd_(i), "(?<sub>OAS3\d{4})_MR_(?<ses>d\d{4})", "names");
                if any(contains(mgb, re.sub))
                    mgd(j) = mgd_(i); %#ok<AGROW>
                    j = j + 1;
                end
            end            
            popd(pwd0); 
        end
        function mgfm = mglob_fieldmap(this)
            mgb = this.mglob_bold();
            pwd0 = pushd("/data/oasis/OASIS3/MR/fieldmap");
            mgfm_ = mglob("OAS3*_MR_d*");
            mgfm_ = strip(mgfm_, filesep);
            j = 1;
            for i = 1:length(mgfm_)
                re = regexp(mgfm_(i), "(?<sub>OAS3\d{4})_MR_(?<ses>d\d{4})", "names");
                if any(contains(mgb, re.sub))
                    mgfm(j) = mgfm_(i); %#ok<AGROW>
                    j = j + 1;
                end
            end            
            popd(pwd0);    
        end
        function mgf = mglob_FLAIR(this)
            mgb = this.mglob_bold();
            pwd0 = pushd("/data/oasis/OASIS3/MR/FLAIR");
            mgt1_ = mglob("OAS3*_MR_d*");
            mgt1_ = strip(mgt1_, filesep);
            j = 1;
            for i = 1:length(mgt1_)
                re = regexp(mgt1_(i), "(?<sub>OAS3\d{4})_MR_(?<ses>d\d{4})", "names");
                if any(contains(mgb, re.sub))
                    mgf(j) = mgt1_(i); %#ok<AGROW>
                    j = j + 1;
                end
            end            
            popd(pwd0);    
        end
        function mgs = mglob_surf(this)
            mgb = this.mglob_bold();
            pwd0 = pushd("/data/oasis/OASIS3/FreeSurfer");
            mgs_ = mglob("OAS3*_Freesurfer53_d*");
            mgs_ = strip(mgs_, filesep);
            j = 1;
            for i = 1:length(mgs_)
                re = regexp(mgs_(i), "(?<sub>OAS3\d{4})_Freesurfer53_(?<ses>d\d{4})", "names");
                if any(contains(mgb, re.sub))
                    mgs(j) = mgs_(i); %#ok<AGROW>
                    j = j + 1;
                end
            end            
            popd(pwd0);    
        end
        function mgt1 = mglob_T1w(this)
            mgb = this.mglob_bold();
            pwd0 = pushd("/data/oasis/OASIS3/MR/T1w");
            mgt1_ = mglob("OAS3*_MR_d*");
            mgt1_ = strip(mgt1_, filesep);
            j = 1;
            for i = 1:length(mgt1_)
                re = regexp(mgt1_(i), "(?<sub>OAS3\d{4})_MR_(?<ses>d\d{4})", "names");
                if any(contains(mgb, re.sub))
                    mgt1(j) = mgt1_(i); %#ok<AGROW>
                    j = j + 1;
                end
            end            
            popd(pwd0);    
        end
        function mgt2 = mglob_T2w(this)
            mgb = this.mglob_bold();
            pwd0 = pushd("/data/oasis/OASIS3/MR/T2w");
            mgt2_ = mglob("OAS3*_MR_d*");
            mgt2_ = strip(mgt2_, filesep);
            j = 1;
            for i = 1:length(mgt2_)
                re = regexp(mgt2_(i), "(?<sub>OAS3\d{4})_MR_(?<ses>d\d{4})", "names");
                if any(contains(mgb, re.sub))
                    mgt2(j) = mgt2_(i); %#ok<AGROW>
                    j = j + 1;
                end
            end            
            popd(pwd0);    
        end

        function [s,r] = migrate_bold2chpc(this)
            pwd0 = pushd("/data/oasis/OASIS3/MR/bold");
            mg = this.mglob_bold();
            for i = 1:length(mg)
                [s,r] = mysystem(sprintf( ...
                    "rsync -ra %s %s/bold", mg(i), this.sourcedata));
            end       
            popd(pwd0);
        end
        function [s,r] = migrate_dwi2chpc(this)
            pwd0 = pushd("/data/oasis/OASIS3/MR/dwi");
            mgd = this.mglob_dwi();
            for i = 1:length(mgd)
                [s,r] = mysystem(sprintf( ...
                    "rsync -ra %s %s/dwi", mgd(i), this.sourcedata));
            end
            popd(pwd0);
        end
        function [s,r] = migrate_fieldmap2chpc(this)
            pwd0 = pushd("/data/oasis/OASIS3/MR/fieldmap");
            mg = this.mglob_fieldmap();
            for i = 1:length(mg)
                [s,r] = mysystem(sprintf( ...
                    "rsync -ra %s %s/fieldmap", mg(i), this.sourcedata));
            end
            popd(pwd0);
        end
        function [s,r] = migrate_FLAIR2chpc(this)
            pwd0 = pushd("/data/oasis/OASIS3/MR/FLAIR");
            mg = this.mglob_FLAIR();
            for i = 1:length(mg)
                [s,r] = mysystem(sprintf( ...
                    "rsync -ra %s %s/FLAIR", mg(i), this.sourcedata));
            end
            popd(pwd0);
        end
        function [s,r] = migrate_surf2chpc(this)
            pwd0 = pushd("/data/oasis/OASIS3/FreeSurfer");
            mg = this.mglob_surf();
            for i = 1:length(mg)
                [s,r] = mysystem(sprintf( ...
                    "rsync -ra %s %s/FreeSurfer", mg(i), this.sourcedata));
            end
            popd(pwd0);
        end
        function [s,r] = migrate_T1w_to_chpc(this)
            pwd0 = pushd("/data/oasis/OASIS3/MR/T1w");
            mg = this.mglob_T1w();
            for i = 1:length(mg)
                [s,r] = mysystem(sprintf( ...
                    "rsync -ra %s %s/T1w", mg(i), this.sourcedata));
            end
            popd(pwd0);
        end
        function [s,r] = migrate_T2w_to_chpc(this)
            pwd0 = pushd("/data/oasis/OASIS3/MR/T2w");
            mg = this.mglob_T2w();
            for i = 1:length(mg)
                [s,r] = mysystem(sprintf( ...
                    "rsync -ra %s %s/T2w", mg(i), this.sourcedata));
            end
            popd(pwd0);
        end

        %% having migrated OASIS3 data to CHPC, arrange HCP filesystem conventions

        function [s,r] = arrange_oasis2hcp(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
