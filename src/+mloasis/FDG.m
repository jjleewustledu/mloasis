classdef FDG < handle & mladni.FDG
    %% FDG uses object-oriented mlfsl.Flirt & 4dfp for registration and mlfourd.ImagingContext2 for actions on imaging data.  
    %  Most common use case:
    %  >> an_fdg = mlfourd.ImagingContext2(fq_filename)
    %  >> obj = mloasis.FDG(an_fdg);
    %  >> obj.call()
    %  >> obj.view()
    %
    %  sourcedata:  1372 subjects; 2680 total pet, all tracers; 127 FDG pet.
    %  derivatives  111 subjects; 127 FDG pet. 
    %  
    %  Created 18-Oct-2022 14:39:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mloasis/src/+mloasis.
    %  Developed on Matlab 9.13.0.2049777 (R2022b) for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function propcluster_tiny()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '6000'; % in MB; deepmrseg requires 10 GB; else 6 GB
            c.AdditionalProperties.Node = '';
            %c.AdditionalProperties.Partition = 'test';
            c.AdditionalProperties.WallTime = '4:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function propcluster()
            c = parcluster;
            c.AdditionalProperties.EmailAddress = '';
            c.AdditionalProperties.EnableDebug = 1;
            c.AdditionalProperties.GpusPerNode = 0;
            c.AdditionalProperties.MemUsage = '10000'; % in MB; deepmrseg requires 10 GB; else 6 GB
            c.AdditionalProperties.Node = '';
            c.AdditionalProperties.Partition = '';
            c.AdditionalProperties.WallTime = '24:00:00'; % 24 h
            c.saveProfile
            disp(c.AdditionalProperties)
        end
        function [j,c] = parcluster_tiny()
            %% #PARCLUSTER_TINY
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            %  --------------------------------------------------------------
            
            c = parcluster;
            disp(c.AdditionalProperties)

            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'OASIS3', 'bids', 'derivatives', 'globbed.mat');
            if ~isfile(globbing_file)
                warning('mloasis:RuntimeWarning', ...
                    'FDG.parcluster_tiny: globbed.mat needed on client side of parallel server')
                j = c.batch(@mloasis.FDG.batch_globbed, 1, {}, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);
                mlbash(sprintf('rsync -a login3.chpc.wustl.edu:%s %s', ...
                    fullfile('"/scratch/jjlee/Singularity/OASIS3/bids/derivatives/globbed.*"'), ...
                    myfileparts(globbing_file)));
                return
            end     

            ld = load(globbing_file);
            globbed = ld.globbed{1};
            globbed = globbed(1:3);
            j = c.batch(@mloasis.FDG.batch, 1, {'globbed', globbed}, 'Pool', 3, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);      
        end
        function [j,c] = parcluster()
            %% #PARCLUSTER
            %  See also https://sites.wustl.edu/chpc/resources/software/matlab_parallel_server/
            %  --------------------------------------------------------------
            %  Use for call_resolved(), batch_revisit_pve1()
            
            c = parcluster;
            disp(c.AdditionalProperties)
            
            globbing_file = fullfile(getenv('SINGULARITY_HOME'), 'OASIS3', 'bids', 'derivatives', 'globbed.mat');
            disp(globbing_file)
            if ~isfile(globbing_file)
                warning('mloasis:RuntimeWarning', ...
                    'FDG.parcluster_tiny: globbed.mat needed on client side of parallel server')
                j = c.batch(@mloasis.FDG.batch_globbed, 1, {}, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);
                mlbash(sprintf('rsync -a login3.chpc.wustl.edu:%s %s', globbing_file, myfileparts(globbing_file)));
                return
            end            
            ld = load(globbing_file);
            for ji = 1:length(ld.globbed)
                j{ji} = c.batch(@mloasis.FDG.batch, 1, {'globbed', ld.globbed{ji}}, 'Pool', 13, ...
                    'CurrentFolder', '.', 'AutoAddClientPath', false);  %#ok<AGROW>
            end
        end
        function globbed = batch_globbed(varargin)
            %% creates globbed.csv, globbed.mat
            
            ip = inputParser;
            addParameter(ip, 'proc', '*', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            g = globT( ...
                fullfile(getenv('SINGULARITY_HOME'), 'OASIS3', 'bids', 'sourcedata', ...
                    'sub-*', 'ses-*', 'pet', '*FDG_pet.nii.gz'));
            Ncores = 31;
            Njobs = ceil(length(g)/Ncores);
            globbed = cell(1,Njobs);
            for gi = 1:Njobs-1
                globbed{gi} = g((gi-1)*Ncores+1:gi*Ncores);
            end
            globbed{Njobs} = g((Njobs-1)*Ncores+1:end);                
            
            t = cell2table(ascol(g), 'VariableNames', {'rawdata_pet_filename'}); % single col of filenamaes
            writetable( ...
                t, fullfile(getenv('SINGULARITY_HOME'), 'OASIS3', 'bids', 'derivatives', 'globbed.csv'));
            save( ...
                fullfile(getenv('SINGULARITY_HOME'), 'OASIS3', 'bids', 'derivatives', 'globbed.mat'), 'globbed');            
                % cell array of filenames grouped for serial runs in parallel
            fprintf('mloasis.FDG.batch.globbed:\n')
            disp(globbed)
        end
        function t = batch(varargin)
            %% Requires completion of mloasis.batch_globbed().
            %  Params:
            %      len (numeric):  # FDG sessions
            %      proc (text):  default 'CASU'
            
            ip = inputParser;
            addParameter(ip, 'globbed', {}, @iscell)
            addParameter(ip, 'proc', '*', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            disp('Start mloasis.FDG.batch()')        

            globbed = ipr.globbed;          
            
            t0 = tic;
            mladni.CHPC3.clean_tempdir();
            for idx = 1:length(globbed)            
                mladni.CHPC3.setenvs();
                try
                    setenv('DEBUG', '')
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    this = mloasis.FDG(fdg_);
                    if ~isempty(this.t1w)
                        this.call_resolve(); 
                        %this.finalize(); 
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end
            t = toc(t0);

            disp('mloasis.FDG.batch() completed')            
        end
        function t = run_local(globbed)
            %% Run interactively on server.  Not for parcluster.
            %  Args:

            arguments
                globbed {mustBeText}
            end
            
            t0 = tic;
            parfor (idx = 1:length(globbed), 16)  
                try
                    setenv('DEBUG', '')
                    fdg_ = mlfourd.ImagingContext2(globbed{idx});
                    this = mloasis.FDG(fdg_);
                    if ~isempty(this.t1w)
                        %this.call_resolve(); 
                        this.build_fdg_renormalized();
                    end
                catch ME
                    handwarning(ME)
                    disp(struct2str(ME.stack))
                end
            end
            t = toc(t0);          
        end
        function missing = find_missing(to_find, locations)
            arguments                
                to_find {mustBeTextScalar} = "*_dlicv.nii.gz"
                locations {mustBeText} = mglob("sub-*/ses-*/pet")
            end

            assert(strcmp(mybasename(pwd), "derivatives"))

            missing = strings(size(locations));
            for i = 1:length(locations)
                mg = mglob(fullfile(locations(i), to_find));
                if isempty(mg)
                    missing(i) = locations(i); end
            end
            missing(missing=="") = [];
        end
        function ic = timeAveragedAndBlurred(ic)
            ic_ = copy(ic);
            t = readtable(strcat(ic.fqfileprefix, '.tsv'), 'FileType', 'text');
            taus = t.Frame_duration;
            taus(1:end-6) = 0;
            ic = ic.timeAveraged('taus', taus);
            ic = ic.blurred(5); % blur ECAT EXACT HR+ to reach 8 mm fwhm standard of ADNI
            ic.filepath = strrep(ic.filepath, 'sourcedata', 'derivatives');
            ensuredir(ic.filepath)
            save(ic);
            mladni.FDG.jsonrecode(ic_, [], ic);
        end
    end

    methods
        % function this = build_fdg_warped(this)
        %     fdgw = mlfourd.ImagingContext2( ...
        %         this.ants.antsApplyTransforms(this.fdg_on_t1w, this.atl_brain, this.t1w_brain));
        %     fdgw = fdgw.thresh(0);
        %     this.fdg_detJ_warped_ = fdgw;
        %     this.fdg_detJ_warped_.fqfp = strcat(this.fdg_on_t1w.fqfp, '_nodetJ_Warped');
        %     this.fdg_detJ_warped_.save();
        % 
        %     mladni.FDG.jsonrecode(...
        %         fdgw, ...
        %         struct('state_changes', 'fdgw = fdgw.thresh(0); this.fdg_detJ_warped_ = fdgw', ...
        %                'image_activity', this.image_mass(this.fdg_detJ_warped_)), ...
        %         this.fdg_detJ_warped_);
        % 
        %     deleteExisting(fdgw);
        % end
        function this = build_fdg_renormalized(this)
            %% renormalizes all FDG by {pons, cerebellar_vermis} per Susan Landau's methods
            
            try
                ic_ = this.fdg_warped_dlicv;
                ic_.selectNiftiTool();
                fp = ic_.fileprefix;
                S = struct('suvr_pons_vermis', nan, 'icv', nan);

                if contains(this.fdg_reference, 'ponsvermis')
                    suvr = this.suvr_pons_vermis();
                    %assert(suvr > 0.5 && suvr < 2);
                    ic_ = ic_./suvr;
                    S.suvr_pons_vermis = suvr;
                end
                if contains(this.fdg_reference, 'icv')
                    icv = this.icv_from_imaging(this.t1w_dlicv);                    
                    ic_ = ic_./icv;
                    S.icv = icv;
                end

                ic_.fileprefix = strcat(fp, '_', this.fdg_reference);
                ic_.save();  
                mladni.FDG.jsonrecode(this.fdg_warped_, S, ic_);
                this.fdg_warped_dlicv_ = ic_;
            catch ME % suvr not always available
                handwarning(ME)
            end
        end
        function this = call_recon_all(this)
            re = regexp(this.t1w.filepath, '(?<subs_dir>\S+)/(?<sub>sub-OAS3\d+)/(?<ses>ses-d\d+)', 'names');
            subs_dir = re.subs_dir;
            sub = re.sub;
            setenv('SUBJECTS_DIR', fullfile(subs_dir, sub));
            mlbash(sprintf('recon-all -s %s -i %s', ses, this.t1w.fqfn));
        end
        function this = call_resolve(this)
            %% #CALL_RESOLVE
            
            %this = this.N4BiasFieldCorrection();
            %this = this.deepmrseg_apply_docker();
            this = this.CreateJacobianDeterminantImages();
            this = this.build_deepmrseg_warped();
            
            this = this.fast();
            this = this.build_fast_warped2();

            this = this.resolve_fdg2t1w();
            if max(this.resolve_err(this.fdg_on_t1w)) > 8
                this = this.flirt_fdg2t1w();
            end
            this = this.build_fdg_warped();
            this.fdg_warped_ = this.repackage_in_atl(this.fdg_warped);
            this = this.build_fdg_warped_masked();
            this.fdg_warped_dlicv_ = this.repackage_in_atl(this.fdg_warped_dlicv);
            this = this.build_fdg_renormalized();
            %this = this.save_qc();
        end


        function this = deepmrseg_apply_docker(this)
            %% e.g.:
            %  nvidia-docker run -it -v $(pwd):/scratch --rm jjleewustledu/deepmrseg_image:20220615 --task dlicv --inImg bert_T1.nii.gz --outImg bert_T1_DLICV.nii.gz

            if isfile(this.t1w_dlicv.fqfn)
                return
            end
            
            %% dlicv

            cmd = sprintf("nvidia-docker run -it -v %s:/scratch --rm jjleewustledu/deepmrseg_image:20220615 --task dlicv --inImg %s --outImg %s", ...
                this.t1w.filepath, this.t1w.filename, this.t1w_dlicv.filename);
            mlbash(cmd);
            mladni.FDG.jsonrecode( ...
                this.t1w, ...
                struct('bash', cmd, ...
                       'image_mass', this.image_mass(this.t1w_dlicv)), ...
                this.t1w_dlicv);
        end
        function ic1 = repackage_in_atl(this, ic)
            ifc = this.atl.imagingFormat;
            ifc.img = ic.imagingFormat.img;
            ifc.fqfp = ic.fqfp;
            ic1 = mlfourd.ImagingContext2(ifc);
            ic1.save();
        end
        function this = resolve_fdg2t1w(this)
            pwd0 = pushd(this.fdg.filepath);
            
            msks{1} = this.t1w_dlicv;
            msks{2} = mlfourd.ImagingContext2('none.nii.gz');
            imgs{1} = this.t1w_blurred;
            imgs{2} = this.fdg;
            t4rb = mlfourd.SimpleT4ResolveBuilder( ...
                'workpath', this.fdg.filepath, ...
                'maskForImages', msks, ...
                'theImages', imgs, ...
                'debug', this.debug);
            t4rb = t4rb.resolve();
            
            movefile(this.niigz(t4rb.theImagesFinal{2}), this.niigz(this.fdg_on_t1w));
            movefile(this.json(t4rb.theImagesFinal{2}), this.json(this.fdg_on_t1w));    
            mladni.FDG.jsonrecode( ...
                this.fdg_on_t1w, ...
                struct('image_activity', this.image_mass(this.fdg_on_t1w)), ...
                this.fdg_on_t1w);
            
            if ~this.debug
                t4rb.deleteFourdfp(t4rb.theImages);
                t4rb.deleteFourdfp(t4rb.theImagesOp(:,1));
            end

            popd(pwd0);
        end
        function s    = suvr_pons_vermis(this)
            pet = this.fdg_warped;
            msk = mlfourd.ImagingContext2( ...
                fullfile(getenv("REFDIR"), "ponsvermis_on_MNI152_T1_2mm.nii.gz"));
            s = pet.volumeAveraged(msk);
            s = s.imagingFormat.img;
        end
        function this = FDG(fdg)
            %% FDG 
            %  Args:
            %      fdg (required any): understood by mlfourd.ImagingContext2.  rawdata gets copied to derivatives.
            %      t1w (any): understood by mlfourd.ImagingContext2.
            %      atl (any): param understood by mlfourd.ImagingContext2.
            %      rawdata_folder(text)
            %      rawdata_path (text): path to replace path containing sub-*/ses-*/{pet,anat}.  Default := fdg.
            %      derivatives_path (text): path to replace path containing sub-*/ses-*/pet.  Default := fdg.
            %      blur (scalar): default := 0.
            %      proc (text): default := 'CASU'.
            %
            %      construct blur
            %      https://doi.org/10.1016/j.nima.2006.08.110 for specs of EE HR+
            %      sqrt(5.65^2 + 4.64^2 + 5.33^2) = 9.047706891804133 % 10 cm radial position
            %      sqrt(4.82^2 + 4.39^2 + 5.1^2) = 8.277348609307209 % 1 cm radial
            %
            %      Bob Koeppe's table "ADNI_PET_ScannerSmoothingTable.xls however, lists 5 mm FWHM isotropic smoothing
            %      to take ECAT EXACT HR+ 6 mm to ADNI standard of 8 mm.
            
            this = this@mladni.FDG(fdg, fdg_proc='CASU');
            
            % FDGDemographics
            this.demog_ = mloasis.FDGDemographics();
            this.fdg_reference_ = 'ponsvermis';

            if ~isempty(getenv('DEBUG'))
                disp(this)
            end
        end
    end

    %% PROTECTED

    methods (Access = protected)
        function construct_paths(this, ipr)
            if 0 == strlength(ipr.rawdata_path)
                re = regexp(this.fdg_.filepath, '(?<rdp>\S+)/sub-OAS3\d+/ses-d\d+/pet', 'names');
                this.rawdata_path_ = re.rdp;
                this.rawdata_path_ = strrep(this.rawdata_path_, 'derivatives', ipr.rawdata_folder);
                assert(isfolder(this.rawdata_path_))
            end
            if 0 == strlength(ipr.derivatives_path)
                this.derivatives_path_ = fullfile(fileparts(this.rawdata_path_), 'derivatives', '');
                assert(~isempty(this.derivatives_path_))
            end
        end
        function fqfn = findT1w(this, ~)
            fqfn = '';

            % check if T1w already found and placed alongside fdg
            deleteExisting(fullfile(this.fdg.filepath, '*_T1w_reslice*'));
            g = glob(fullfile(this.fdg.filepath, '*_T1w.nii.gz'));
            if ~isempty(g)
                fqfn = g{1};
                return
            end
            
            try
                % find /pth/to/rawdata/sub-123S4566/ses-yyyymmdd/anat/*-T1w.nii.gz
                fdgpth_ = strrep(this.fdg.filepath, 'derivatives', 'sourcedata'); % /pth/to/sourcedata/sub-OAS30001/ses-d1234/pet
                sesfold_ = mybasename(myfileparts(fdgpth_)); % ses-d1234
                re = regexp(sesfold_, 'ses-d(?<dt>\d+)', 'names');
                dt_fdg_ = str2double(re.dt);

                subpth_ = fileparts(myfileparts(fdgpth_)); % /pth/to/rawdata/sub-123S4566
                globt1w_ = globT(fullfile(subpth_, 'ses-*', 'anat', '*_T1w.nii.gz'));
                globt1w_ = globt1w_(~contains(globt1w_, 'mask'));
                %if any(contains(globt1w_, 'irspgr')) % KLUDGE for irspgr containing t1w and t1w containing localizer                    
                %    globt1w_ = globt1w_(contains(globt1w_, 'irspgr'));
                %end
                if isempty(globt1w_)
                    return
                end                
                
                % find globt1w__ for T1w datetime closest to FDG datetime
                dts_ = nan(size(globt1w_));
                for ig = 1:length(globt1w_)
                    sesfold_ = mybasename(myfileparts(myfileparts(globt1w_{ig})));
                    re = regexp(sesfold_, 'ses-d(?<dt>\d+)', 'names');
                    dts_(ig) = str2double(re.dt);
                end
                if isempty(dts_)
                    return
                end
                ddt = min(abs(dts_ - dt_fdg_)); % ddt of anat nearest in time to fdg
                select = abs(dts_ - dt_fdg_) == ddt;
                globt1w__ = globt1w_(select);
                
                % find idx_proc of proc_ with longest processing descriptor
                proc_ = cell(size(globt1w__));
                for ip = 1:length(globt1w__)
                    [~,fp_] = myfileparts(globt1w__{ip});
                    re = regexp(fp_, 'sub-OAS3\d+_ses-d\d+(?<proc>\S+)T1w', 'names');
                    proc_{ip} = re.proc;
                end
                [~,idx_proc] = max(cell2mat(cellfun(@length , proc_, 'UniformOutput', false))); % idx of longest proc specification
                if isempty(idx_proc)
                    return
                end

                % =======================================
                if ~isempty(getenv('DEBUG'))
                    fdgpth_ %#ok<*NOPRT> 
                    sesfold_
                    dt_fdg_
                    subpth_
                    fprintf('mladni.FDG.findT1w.globpth_:\n')
                    dts_'
                    fprintf('mladni.FDG.findT1w.globt1w_:\n');
                    disp(ascol(globt1w_))
                    fprintf('mladni.FDG.findT1w.proc_:\n');
                    disp(ascol(proc_))
                    idx_proc
                end
                % =======================================
                fqfn = globt1w__{idx_proc};
                if ~isfile(fqfn)
                    return
                end
                if ~isempty(getenv('DEBUG'))
                    fprintf('mladn.FDG.findT1w:  found %s\n', fqfn);
                end
            catch ME
                handwarning(ME)
            end
        end
        function icd  = prepare_derivatives(this, ic)

            if isempty(this.fdg_)
                % case fdg, etc.
                deriv_pth = strrep(ic.filepath, 'rawdata', 'derivatives');
            else
                % case t1w
                deriv_pth = this.fdg_.filepath;
            end
            ensuredir(deriv_pth);
            try
                if ~isfile(fullfile(deriv_pth, ic.filename))
                    mysystem(sprintf('cp -f %s %s', ic.fqfn, deriv_pth), '-echo');
                end
                if ~isfile(fullfile(deriv_pth, strcat(ic.fileprefix, '.json'))) && ...
                        isfile(strcat(ic.fqfp, '.json'))
                    mysystem(sprintf('cp -f %s %s', strcat(ic.fqfp, '.json'), deriv_pth), '-echo');
                end
                if ~isfile(fullfile(deriv_pth, strcat(ic.fileprefix, '.tsv'))) && ...
                        isfile(strcat(ic.fqfp, '.tsv'))
                    mysystem(sprintf('cp -f %s %s', strcat(ic.fqfp, '.tsv'), deriv_pth), '-echo');
                end
            catch
            end
            icd = mlfourd.ImagingContext2(fullfile(deriv_pth, ic.filename));

            % OASIS FDG is dynamic, so create static
            if contains(icd.fileprefix, "FDG") && ...
                    ~contains(icd.fileprefix, "_avgt_b50")
                icd = this.timeAveragedAndBlurred(icd);
            end
            if ~isfile(icd)
                icd.save();
            end
            % this.fdg_ = icd; % bug that assigns this.fdg_ <- t1w

            mlbash(sprintf('chmod -R 755 %s', deriv_pth));
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
