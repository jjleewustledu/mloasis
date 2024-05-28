classdef OasisDemographics < handle
    %% 
    % __Tom Earnest 12:00__
    %
    % Relevant post I had a bit ago about amyloid/clinical measures having missing data in some places but not others:
    % https://cil-wustl.slack.com/archives/C02NDGL984T/p1664574205433739
    % 
    % __Tom Earnest__ Last post of the day :woozy_face: 
    % 
    % For those working with OASIS3 data, I have ran into a couple instances
    % where it looked like there was a fair amount of missing data in one table, but elsewhere I found a table with
    % seemingly complete values.  I'm not quite sure what the discrepancy is (e.g. one table hasn't been updated / one table
    % is correct / one table is filtered for QC or otherwise / the two tables actually contain different information).  And
    % I might be missing something obvious - but anyway here are the cases I encountered in case this helps anyone:
    % 
    % __Centiloids / amyloid positivity:__ On XNAT central, if you go to Browse > My Projects > OASIS3, then Add Tab > PUPs, you
    % will get a table that has centiloids computed for binding potential/SUVR & with/without PVC correction.  All 4
    % measures have substantial missing values.  But if you do Options > Edit Columns, you can add many more fields, one
    % being "PET_fSUVR_rsf_TOT_CORTMEAN" (near the bottom).  This column has complete values, and can be converted to
    % centiloids using the equations in the OASIS3 documentation.  Computing centiloids from this column gave a value that
    % was very, very near the centiloids that are shown by default (though not precisely the same).  I can share a figure of
    % this if anyone is curious. 
    % 
    % __CDR / Demographics:__ If you go to Browse > My Projects > OASIS3, then Add Tab > ADRC
    % Clinical Data you will get a table which has a "cdr" column.  This appears to be complete, but does not have entries
    % for all subjects (e.g. OASIS30055, one of the subjects with tau imaging).   But if you go to Browse > My Projects >
    % OASIS3 and click the subject named "0AS_data_files" > MR Session (OASIS3_data_files), you get to a page where you can
    % download several CSV files.  The "UDSb4" (Form B4: Global Staging CDR: Standard and Supplemental) has CDR values which
    % cover the the tau subset (as well as everyone else, apparently).  There is also "demo" which has complete demographic
    % information. 
    % 
    % Show less Thread in minds_ad | Sep 30th, 2022 | View message
    %  
    %  Created 20-Apr-2023 21:35:49 by jjlee in repository /Users/jjlee/MATLAB-Drive/mloasis/src/+mloasis.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        cutoff_pib
        cutoff_av45
        days_separation_tol
        workdir
    end

    methods
        function this = OasisDemographics(varargin)
            this.workdir = fullfile(getenv("MATLABDRIVE"), "mloasis", "data");
            this.cutoff_pib = 16.4; % CL PIB MCSURVR RSF, OASIS-3_Imaging_Data_Dictionary_v2.3.pdf
            this.cutoff_av45 = 20.6; % CL AV45 MCSURVR RSF OASIS-3_Imaging_Data_Dictionary_v2.3.pdf
            this.days_separation_tol = 365;
        end

        function T = age_info(this)
            warning("off", "MATLAB:table:ModifiedAndSavedVarnames")
            T_ = readtable(fullfile(this.workdir, "OASIS3_UDSb4_cdr.csv"));
            sub = cell2mat(cellfun(@get_sub, T_.OASISID, UniformOutput=false));
            ses = T_.days_to_visit;
            age = T_.ageAtVisit;
            T = natsortrows(table(sub, ses, age), [], [1, 2]);
            
            T(T.ses < -36500, :) = []; % delete mangled data
            warning("on", "MATLAB:table:ModifiedAndSavedVarnames")

            function s = get_sub(s_)
                re = regexp(s_, "OAS(?<sub>\d{5})$", "names");
                s = str2double(re.sub);
            end
        end
        
        function T = apoe4_info(this)
            warning("off", "MATLAB:table:ModifiedAndSavedVarnames")
            T_ = readtable(fullfile(this.workdir, "oasis3_clinical_full.csv"));
            sub = cell2mat(cellfun(@get_sub, T_.ADRC_ADRCCLINICALDATAID, UniformOutput=false));

            unique_subs = ascol(unique(sub));
            apoe4 = NaN(size(unique_subs));
            for u = 1:length(unique_subs)
                select_sub = sub == unique_subs(u);
                U_ = T_(select_sub, :);
                apoe4(u) = count(num2str(U_.apoe(1)), '4');
            end
            T = natsortrows(table(unique_subs, apoe4));
            warning("on", "MATLAB:table:ModifiedAndSavedVarnames")

            function s = get_sub(s_)
                try
                    re = regexp(s_, "OAS(?<sub>\d{5})_ClinicalData_d(?<ses>\d+)$", "names");
                    s = str2double(re.sub);
                catch ME
                    handwarning(ME)
                    s = NaN;
                end
            end
        end
        
        function T = sex_info(this) 
            warning("off", "MATLAB:table:ModifiedAndSavedVarnames")
            T_ = readtable(fullfile(this.workdir, "oasis3_clinical_full.csv"));
            sub = cell2mat(cellfun(@get_sub, T_.ADRC_ADRCCLINICALDATAID, UniformOutput=false));
            
            unique_subs = ascol(unique(sub));
            sex = cell(size(unique_subs));
            for u = 1:length(unique_subs)
                select_sub = sub == unique_subs(u);
                U_ = T_(select_sub, :);
                sex{u} = U_.M_F{1};
            end
            sex = ascol(categorical(sex));
            T = natsortrows(table(unique_subs, sex));            
            warning("on", "MATLAB:table:ModifiedAndSavedVarnames")

            function s = get_sub(s_)
                try
                    re = regexp(s_, "OAS(?<sub>\d{5})_ClinicalData_d(?<ses>\d+)$", "names");
                    s = str2double(re.sub);
                catch ME
                    handwarning(ME)
                    s = NaN;
                end
            end
        end

        function T = amyloid_info(this)
            T_ = this.centiloid_info();
            sub = T_.sub;
            ses = T_.ses;

            amyloidosis = false(size(sub)); % conceptually easier to "opt in" amyloidosis by rules of get_amyloidosis()
            unique_subs = asrow(unique(sub));
            for asub = unique_subs
                select_sub = sub == asub;
                U_ = T_(select_sub, :);
                amy_ = get_amyloidosis(U_);
                amyloidosis(select_sub) = amy_;
            end

            T = natsortrows(table(sub, ses, amyloidosis), [], [1,2]);

            function amy = get_amyloidosis(U__)
                amy = ...
                    (strcmp(U__.tracer, 'AV45') & U__.centiloid > this.cutoff_av45) | ...
                    (strcmp(U__.tracer, 'PIB') & U__.centiloid > this.cutoff_pib);
            end
        end
        
        function call(this)
            %% prepare for running NMF

            assert(~isempty(this.table_fdg));
            t = this.table_cn(true);
            Filelist = strrep(t.Filelist, '/home/usr', '/scratch');
            u = table(Filelist);
            fqfn = fullfile(this.workdir, 'baseline_cn', 'nifti_files.csv');
            ensuredir(myfileparts(fqfn));
            writetable(u, fqfn, WriteVariableNames=false);
        end
        
        function T = cdr_info(this)
            warning("off", "MATLAB:table:ModifiedAndSavedVarnames")
            T_ = readtable(fullfile(this.workdir, "OASIS3_UDSb4_cdr.csv"));
            sub = cell2mat(cellfun(@get_sub, T_.OASISID, UniformOutput=false));
            ses = T_.days_to_visit;
            cdr = T_.CDRTOT;
            T = natsortrows(table(sub, ses, cdr), [], [1, 2]);
            
            T(T.ses < -36500, :) = []; % delete mangled data
            warning("on", "MATLAB:table:ModifiedAndSavedVarnames")

            function s = get_sub(s_)
                re = regexp(s_, "OAS(?<sub>\d{5})$", "names");
                s = str2double(re.sub);
            end
        end
        
        function T = centiloid_info(this)
            T_ = readtable(fullfile(this.workdir, "OASIS3_amyloid_centiloid.csv"));
            sub = cellfun(@(x) str2double(x(4:end)), T_.subject_id);
            ses = cellfun(@(x) str2double(get_ses_days(x)), T_.oasis_session_id);
            tracer = T_.tracer;
            centiloid = T_.Centiloid_fSUVR_rsf_TOT_CORTMEAN;

            T = natsortrows(table(sub, ses, tracer, centiloid), [], [1, 2]);

            function s = get_ses_days(s_)
                re = regexp(s_, "OAS\d{5}_[A-Z0-9]+_d(?<days>\d{4})", "names");
                s = re.days;
            end
        end
        
        function T = globbed_info(this)
            ld = load(fullfile(this.workdir, "globbed126.mat"));
            globbed = ascol(ld.globbed126);
            sub = NaN(length(globbed), 1);
            ses = NaN(length(globbed), 1);
            for idx = 1:length(globbed)
                re = regexp(mybasename(globbed(idx)), "sub-OAS(?<sub>\d{5})_ses-d(?<ses>\d{4})_\S+", "names");
                sub(idx) = str2double(re.sub);
                ses(idx) = str2double(re.ses);
            end
            Filelist = myfileprefix(globbed) + "_avgt_b50_on_T1w_Warped_dlicv_ponsvermis.nii.gz";
            T = natsortrows(table(Filelist, globbed, sub, ses), [], [2,3]);
        end
        
        function T = table_fdg(this)
            %% globbed_info with matched amyloid, cdr, and considering days_separation
            %  N = 108

            T_ = this.globbed_info;
            globbed = T_.globbed;
            Filelist = T_.Filelist;
            sub = T_.sub;
            ses = T_.ses;
            Tcdr = this.cdr_info();
            Tamy = this.amyloid_info();
            Tage = this.age_info();
            Tsex = this.sex_info();
            Tapoe4 = this.apoe4_info();

            warning("off", "MATLAB:badsubscript")
            cdr = NaN(size(sub));
            amyloidosis = NaN(size(sub));
            age = NaN(size(sub));
            sex = cell(size(sub));
            apoe4 = NaN(size(sub));
            for row = 1:size(T_,1)
                the_sub = sub(row);
                the_ses = ses(row);

                try % cdr
                    Ucdr = Tcdr(Tcdr.sub == the_sub, :);
                    dday = abs(Ucdr.ses - the_ses);
                    Ucdr = addvars(Ucdr, dday, NewVariableNames="dday");
                    Ucdr = sortrows(Ucdr, "dday");
                    if Ucdr.dday(1) < this.days_separation_tol
                        cdr(row) = Ucdr.cdr(1);
                    else
                        fprintf("sub %g, ses %g, cdr %g had delta days ~ %g\n", ...
                            the_sub, the_ses, Ucdr.cdr(1), Ucdr.dday(1))
                    end
                catch ME
                    handwarning(ME)
                end

                try % age
                    Uage = Tage(Tage.sub == the_sub, :);
                    dday = the_ses - Uage.ses;
                    Uage = addvars(Uage, dday, NewVariableNames="dday");
                    abs_dday = abs(Uage.ses - the_ses);
                    Uage = addvars(Uage, abs_dday, NewVariableNames="abs_dday");
                    Uage = sortrows(Uage, "abs_dday");
                    age(row) = Uage.age(1) + years(days(Uage.dday(1)));
                catch ME
                    handwarning(ME)
                end

                try % sex
                    Usex = Tsex(Tsex.unique_subs == the_sub, :);
                    sex{row} = Usex.sex(1);
                catch ME
                    handwarning(ME)
                end

                try % apoe4
                    Uapoe4 = Tapoe4(Tapoe4.unique_subs == the_sub, :);
                    apoe4(row) = Uapoe4.apoe4(1);
                catch ME
                    handwarning(ME)
                end

                try % amyloid, must be last task block in for-loop
                    Uamy = Tamy(Tamy.sub == the_sub, :);
                    if ~isempty(Uamy)
                        Uamy__ = Uamy; % chronological
                        dday = abs(Uamy.ses - the_ses);
                        Uamy = addvars(Uamy, dday, NewVariableNames="dday");
                        Uamy = sortrows(Uamy, "dday"); % ordered by abs delta day from the_ses
                        if Uamy.dday(1) < this.days_separation_tol % select any amloid_s1 within datetime_separation_tol of acqdate
                            amyloidosis(row) = Uamy.amyloidosis(1);
                            continue
                        end
                        if all(~Uamy__.amyloidosis) % no known amyloidosis through end of amy scanning
                            if the_ses <= max(Uamy__.ses) + this.days_separation_tol
                                amyloidosis(row) = false;
                                continue
                            else
                                continue
                            end
                        end
                        if Uamy__.amyloidosis(1) % known amyloidosis since start of amy scanning
                            if the_ses >= min(Uamy__.ses) - this.days_separation_tol
                                amyloidosis(row) = true;
                                continue
                            else
                                continue
                            end
                        end

                        % at least one amy+ scan after first amy scan
                        [~,idx] = max(Uamy__.amyloidosis);
                        amyloidosis(row) = the_ses >= Uamy__.ses(idx) - this.days_separation_tol;
                    end
                catch ME
                    handwarning(ME)
                end
            end
            warning("on", "MATLAB:badsubscript")

            %sex = categorical(sex);

            T = table(Filelist, globbed, sub, ses, cdr, amyloidosis, age, sex, apoe4); % age, sex, apoe4
            T = T(~isnan(cdr) & ~isnan(amyloidosis), :);
            T.amyloidosis = logical(T.amyloidosis);
        end
        
        function T = table_ad(this, cs)
            %% N = 5 (long), 4 (cs)

            arguments
                this mloasis.OasisDemographics
                cs logical = false
            end

            T = this.table_fdg();
            T = T(T.cdr >= 0.5 & T.amyloidosis, :);

            if cs
                T = this.table_firstscan(T);
            end
        end
        
        function T = table_cn(this, cs)
            %% N = 68 (long), 63 (cs)

            arguments
                this mloasis.OasisDemographics
                cs logical = false
            end

            T = this.table_fdg();
            T = T(T.cdr == 0 & ~T.amyloidosis, :);

            if cs
                T = this.table_firstscan(T);
            end
        end
        
        function T = table_preclinical(this, cs)
            %% N = 30 (long), 26 (cs)

            arguments
                this mloasis.OasisDemographics
                cs logical = false
            end

            T = this.table_fdg();
            T = T(T.cdr == 0 & T.amyloidosis, :);

            if cs
                T = this.table_firstscan(T);
            end
        end
        
        function T = table_cdr_0p5_apos(this, cs)
            %% N = 4, 4

            arguments
                this mloasis.OasisDemographics
                cs logical = false
            end

            T = this.table_fdg();
            T = T(T.cdr == 0.5 & T.amyloidosis, :);

            if cs
                T = this.table_firstscan(T);
            end
        end
        
        function T = table_cdr_gt_0p5_apos(this, cs)
            %% N = 1, 1
            
            arguments
                this mloasis.OasisDemographics
                cs logical = false
            end

            T = this.table_fdg();
            T = T(T.cdr > 0.5 & T.amyloidosis, :);

            if cs
                T = this.table_firstscan(T);
            end
        end
        
        function T = table_firstscan(~, T_)
            select_cols = contains(T_.Properties.VariableNames, 'sub', IgnoreCase=true) | ...
                contains(T_.Properties.VariableNames, 'ses', IgnoreCase=true);
            T_ = natsortrows(T_, [], select_cols);

            unique_subs_ = unique(T_.sub);
            T = T_(T_.sub == unique_subs_(1), :); 
            T = T(T.ses == min(T.ses), :); 
            for tidx = 2:length(unique_subs_)
                U = T_(T_.sub == unique_subs_(tidx), :); % pick subject
                U = U(U.ses == min(U.ses), :); % pick subject's first scan 
                T = [T; U]; %#ok<AGROW> % append 
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
