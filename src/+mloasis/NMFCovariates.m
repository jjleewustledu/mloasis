classdef NMFCovariates < mladni.NMFCovariates & handle
    %% line1
    %  line2
    %  
    %  Created 13-Apr-2024 13:07:43 by jjlee in repository /Users/jjlee/MATLAB-Drive/mloasis/src/+mloasis.
    %  Developed on Matlab 24.1.0.2537033 (R2024a) for MACA64.  Copyright 2024 John J. Lee.

    methods
        function this = NMFCovariates(varargin)
            this = this@mladni.NMFCovariates(varargin{:});

            this.data_home_ = fullfile(getenv("SINGULARITY_HOME"), "OASIS3");
            this.demogr_ = mloasis.OasisDemographics();
            this.pet_on_T1w_suffix_ = 'orient-rpi_pet_avgt_b50_on_T1w.nii.gz';
        end

        function t = AddAcqDuration(~, t)
            %% AcqDuration ~ floating-point years since first scan
            %  Age ~ adjusted to precision of AcqDuration after first scan

            nans = nan(size(t, 1), 1);
            t = addvars(t, nans, NewVariableNames={'AcqDuration'});
            subs = unique(t.sub);
            for s = asrow(subs)
                sub_select = t.sub == s;
                u = t(sub_select, :);  
                u_acqdur = u.ses;
                u_acqdur = u_acqdur - min(u_acqdur);                
                t.AcqDuration(sub_select) = u_acqdur;
            end
        end

        function t = apply_table_qc(this, t)
            %% qc that is more stringent than that for NMF, needed for regressions

            debug_file = fullfile(this.componentDir, stackstr()+this.datestr()+".mat");
            save(debug_file, 't');

            % start with 3478 rows
            
            vns = t.Properties.VariableNames;
            if any(contains(vns, 'cdr'))
                t = t(t.cdr ~= -1, :);
            end % 3458 rows remaining
            if any(contains(vns, 'Cohort'))
                t = t(t.Cohort ~= 'unknown', :);
            end % 1961 rows remaining

            if any(contains(vns, 'Components_1'))
                t = t(~isnan(t.Components_1), :);
            end % 1961 rows remaining
            if any(contains(vns, 'age'))
                t = t(~isnan(t.age), :);
            end % 1961 rows remaining       

            if any(contains(vns, 'Dlicv'))
                t = t(t.Dlicv > 1e6, :);
            end % 1958 rows remaining
            if any(contains(vns, 'PVE1'))
                t = t(t.PVE1 > 0.3e6, :);
            end % 1958 rows remaining
            if any(contains(vns, 'RegErr'))
                t = t(t.RegErr < 2.5, :); % see also mladni.FDG()
            end % 1941 rows remaining

            if any(contains(vns, 'apoe4'))
                t = t(~isnan(t.apoe4), :);
            end %

            if ~isemptytext(this.EXCLUSIONS)
                for fidx = 1:length(this.EXCLUSIONS)
                    select = ~contains(t.Filelist, this.EXCLUSIONS{fidx});
                    t = t(select, :);
                end
            end
        end

        function t = table_covariates(this, opts)
            %% Saves and returns table useful for inferences using RStudio.
            %  Args:
            %      save_1comp (logical): call this.table_covariates_1comp() for all components.

            arguments
                this mloasis.NMFCovariates
                opts.save_1comp logical = false                
            end

            if ~isempty(this.table_covariates_cache_)
                fprintf("%s: using cached in memory\n", stackstr())
                t = this.table_covariates_cache_;
                return
            end
            cache_file = this.covariates_file;
            if isfile(cache_file)
                fprintf("%s: using cached from filesystem\n", stackstr())
                ld = load(cache_file);
                this.table_covariates_cache_ = ld.t;
                t = this.table_covariates_cache_;
                return
            end

            % from OasisDemographics
            t = table_fdg(this); 

            % Cohort ~ categorical
            Cohort = repmat({'unknown'}, size(t,1), 1);
            Cohort(t.cdr == 0 & t.amyloidosis == 0) = {'CDR=0,amy-'};
            Cohort(t.cdr == 0 & t.amyloidosis == 1) = {'CDR=0,amy+'};
            Cohort(t.cdr == 0.5 & t.amyloidosis == 1) = {'CDR=0.5,amy+'};
            Cohort(t.cdr > 0.5 & t.amyloidosis == 1) = {'CDR>0.5,amy+'};
            Cohort(t.cdr > 0 & t.amyloidosis == 0) = {'CDR>0,amy-'};
            Cohort = categorical(Cohort);
            t = addvars(t, Cohort, NewVariableNames={'Cohort'});

            % DLICV ~ 3478 rows, 12 nans
            % t_ = this.table_dlicv;
            % t = addvars(t, t_.Dlicv, NewVariableNames={'Dlicv'});

            % PVE1 ~ 3478 rows, 8 nans
            t_ = this.table_pve1;
            t = addvars(t, t_.PVE1, NewVariableNames={'PVE1'});

            % RegErr ~ 3478 rows, 38 nans
            t_ = this.table_regerr;
            t = addvars(t, t_.RegErr, NewVariableNames={'RegErr'});

            % Components ~ 3470 rows -> 3478 rows, 8 nans
            % implicitly marks empty Filelist for exclusion by apply_table_qc()
            t_ = this.table_selectedComponentWeightedAverageNIFTI;
            t = this.addvars_by_filelist(t, t_, t_.Components, NewVariableNames={'Components'});
            t = splitvars(t, 'Components');

            % apply table qc
            t = this.apply_table_qc(t);

            % sort rows by Subject, then AcqDate
            % t = sortrows(t, ["Subject", "AcqDate"]);

            % AcqDuration ~ floating-point years since first scan, after removing faulty scans
            t = this.AddAcqDuration(t);

            % store cache, & save/write table
            this.table_covariates_cache_ = t;
            save(cache_file, 't');
            writetable(t, strrep(cache_file, ".mat", ".csv"));

            % save separate tables for each component
            if opts.save_1comp
                for idx = 1:this.selectedNumBases
                    t1 = this.table_covariates_1comp(idx);
                    save(strrep(cache_file, ".mat", "_1comp.mat"), 't1');
                    writetable(t1, strrep(cache_file, ".mat", "_1comp.csv"));
                end
            end
        end

        function t = table_covariates_1stscan(this)
            t = this.demogr_.table_firstscan(this.table_covariates());

            cache_file = this.covariates_1stscan_file;
            save(cache_file, 't');
            writetable(t, strrep(cache_file, ".mat", ".csv"));
        end
        function t = table_fdg(this)
            if ~isempty(this.table_fdg_)
                t = this.table_fdg_;
                return
            end
            this.table_fdg_ = this.demogr_.table_fdg();
            t = this.table_fdg_;
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
