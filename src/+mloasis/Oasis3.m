classdef Oasis3 < mladni.DataCuration & handle
    %% line1
    %  line2
    %  
    %  Created 06-Aug-2023 15:38:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mladni/src/+mladni.
    %  Developed on Matlab 9.12.0.2170939 (R2022a) Update 6 for MACI64.  Copyright 2023 John J. Lee.

    properties (Dependent)
        data_home
        nmf_fdg_home
    end

    methods % GET
        function g = get.data_home(~)
            g = fullfile(getenv("SINGULARITY_HOME"), "OASIS3");
        end
        function g = get.nmf_fdg_home(this)
            g = fullfile(this.data_home, "NMF_FDG");
        end
    end

    methods
        function this = Oasis3()
            this.demogr = mloasis.OasisDemographics();
            this.nmf = mladni.NMF( ...
                data_home=this.data_home, ...
                selectedNumBases=this.selectedNumBases, ...
                numBases=2:2:32);
            this.nmfc = mloasis.NMFCovariates();
            this.nmfh = mladni.NMFHierarchies(data_home=this.data_home);
            this.nmfr = mladni.NMFRadar();
        end

        function anticlust_cn_adni2oasis(this)
            %% Prepares folders in NMF_FDG for reproducibility runs of ADNI to compare to OASIS3.
            %  Requires completion of mladni.Adni.anticlust_cn_repeat()

            R_fqfn = fullfile(this.nmf_fdg_home, "anticlust_oasis_vs_adni.R");
            mysystem(sprintf("Rscript %s", R_fqfn))

            mladni.Anticlust.prepare_folders_for_VolBin(data_home=this.data_home)
        end

        function anticlust_cn_repeat(this)
            %% prepares folders in NMF_FDG for reproducibility runs of NMF

            t = this.demogr.table_cn(true);
            t.Filelist = strrep(t.Filelist, "/home/usr", "/scratch");
            csv_fqfn = fullfile(this.nmf_fdg_home, "anticlust_cn_repeat.csv");
            writetable(t, csv_fqfn, WriteVariableNames=true);

            R_fqfn = fullfile(this.nmf_fdg_home, "anticlust_cn_repeat.R");
            % mysystem(sprintf("Rscript %s", R_fqfn));            

            mladni.Anticlust.prepare_folders_for_VolBin(data_home=this.data_home)
        end

        function call(this)
            % mloasis.FDG ~ computes brain mask (dlicv), warps T1w & FDG to MNI, computes PVE, applies brain mask

            % generate table with filelist for NMF
            % call(this.demogr); 
            
            % NMF on cluster
            % $SINGULARITY_HOME/OASIS3/VolBin/submit_*.sh

            % create anticlust filelists with RStudio, then run on cluster
            this.anticlust_cn_repeat();
            % $SINGULARITY_HOME/ADNI/VolBin/submit_anticlust_20230526.sh

            % write cache "X.mat", then calculate reconstruction errors, reproducibility analysis
            call(this.nmf)

            % calculate component-weighted averages
            call2(this.nmf)

            % build argmax maps
            % for b = 2:2:24
            %     nmfc = mloasis.NMFCovariates(selectedNumBases=b); 
            %     nmfc.table_covariates_1stscan(); 
            % end
            % this.nmfh.build_argmax_maps;

            % check completeness
            % fdg = this.demogr.table_fdg;
            % assert(sum(cellfun(@isempty, fdg.Filelist)) == 0, stackstr())

            % report ARIs, Hungarian overlaps between OASIS3 and ADNI
            this.matrices_ARIs();
            this.matrices_overlaps();

            % create images of means, var, std, median, iqr of baseline_cn, longitudinal_*     
            try
                this.nmf.build_stats_imaging(inputDir=fullfile(this.nmf_fdg_home, "baseline_cn"));
            catch ME
                handwarning(ME)
            end
            % this.nmf.build_table_variances(subgroups);
        end

        function [ari_ao,ari_aa,ari_oo] = matrices_ARIs(this)
            lab = "span " + this.selected_spans;
            adni_cn = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "NMF_FDG", "baseline_cn");
            oasis_cn = fullfile(getenv("SINGULARITY_HOME"), "OASIS3", "NMF_FDG", "baseline_cn");
            ari_aa = this.matrix_ARI(adni_cn, adni_cn);
            f2 = figure;
            h2 = heatmap(lab, lab, ari_aa);
            title("ADNI - ADNI")
            saveFigure2(f2, fullfile(oasis_cn, "ari_aa"))
            ari_oo = this.matrix_ARI(oasis_cn, oasis_cn);
            f3 = figure;
            h3 = heatmap(lab, lab, ari_oo);
            title("OASIS3 - OASIS3")
            saveFigure2(f3, fullfile(oasis_cn, "ari_oo"))
            ari_ao = this.matrix_ARI(adni_cn, oasis_cn);
            f1 = figure;
            h1 = heatmap(lab, lab, ari_ao);
            title("ADNI - OASIS3")
            saveFigure2(f1, fullfile(oasis_cn, "ari_ao"))
        end

        function ari = matrix_ARI(this, pth1, pth2)
            spans = this.selected_spans;
            ari = nan(length(spans), length(spans));
            for si1 = 1:length(spans)
                for si2 = 1:length(spans)
                    n1 = mlfourd.ImagingFormatContext2( ...
                        fullfile(pth1, "NumBases"+spans(si1), "OPNMF", "niiImg", "Basis_argmax_brain_mask.nii"));
                    n2 = mlfourd.ImagingFormatContext2( ...
                        fullfile(pth2, "NumBases"+spans(si2), "OPNMF", "niiImg", "Basis_argmax_brain_mask.nii"));
                    p1 = n1.img(n1.img > 0);
                    p2 = n2.img(n2.img > 0);
                    ari(si1, si2) = mladni.AdjRandIndex.rand_index(p1, p2, 'adjusted');
                end
            end
        end
        
        function overlap_ao = matrices_overlaps(this)
            span = this.selected_spans(end);
            lab = "P" + 1:span;
            adni_cn = fullfile(getenv("SINGULARITY_HOME"), "ADNI", "NMF_FDG", "baseline_cn");
            oasis_cn = fullfile(getenv("SINGULARITY_HOME"), "OASIS3", "NMF_FDG", "baseline_cn");
            overlap_ao = this.matrix_overlap(oasis_cn, adni_cn);
            f1 = figure;
            h1 = heatmap(lab, lab, ari_ao);
            xlabel("ADNI")
            ylabel("OASIS3")
            title("Medians of Hungarian Overlaps")
            saveFigure2(f1, fullfile(oasis_cn, "overlap_ao"));
        end

        function ol = matrix_overlap(this, pth1, pth2)
            span = mladni.NMFHierarchies.selected_spans(end);
            ol = nan(span, span);

            results1 = load(fullfile(pth1, "NumBases"+span, "OPNMF", "ResultsExtractBases.mat"));
            results2 = load(fullfile(pth2, "NumBases"+span, "OPNMF", "ResultsExtractBases.mat"));

            for p1 = 1:span  % pattern_1

                B1 = results1.B(:, p1);

                for p2 = 1:span  % pattern_2

                    B2 = results2.B(:, p2);

                    [~, overlap] = mladni.AdjRandIndex.evaluate_pair(B1, B2, do_ari=false);  % overlap ~ Ncomponents x 1
                    ol(p1, p2) = median(overlap);  % scalar
                end
            end
        end

        function T = table_census(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
