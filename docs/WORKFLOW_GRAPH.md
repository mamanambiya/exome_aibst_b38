# Exome Analysis Workflow Architecture - Visual Graph

## Main Workflow Overview

```mermaid
graph TB
    Start([Start Workflow]) --> Phase1[PHASE 1: Core Data Processing]
    Phase1 --> Phase2[PHASE 2: Population Grouping]
    Phase2 --> Phase3[PHASE 3: Statistical Analyses]
    Phase3 --> Phase4[PHASE 4: HDV Detection]
    Phase4 --> Phase5[PHASE 5: Consequence Analysis]
    Phase5 --> End([Complete])
    
    style Phase1 fill:#e1f5ff
    style Phase2 fill:#fff4e1
    style Phase3 fill:#ffe1f5
    style Phase4 fill:#e1ffe1
    style Phase5 fill:#f5e1ff
```

## PHASE 1: Core Data Processing (data_proc)

```mermaid
graph TB
    subgraph "Input"
        VCF[Raw VCF Files<br/>per Dataset]
        REF[Reference Files<br/>dbSNP, ClinVar, COSMIC]
    end
    
    VCF --> GetMap[get_map<br/>Extract chromosomes]
    GetMap --> GenChunks[generate_chunks_vcf<br/>Create 25MB chunks]
    GenChunks --> SplitChunk[split_vcf_to_chunk<br/>Extract regions]
    
    SplitChunk --> QC[dataset_qc_alt<br/>Quality Control<br/>PASS, biallelic SNPs]
    
    QC --> SnpEff[annotate_snpeff<br/>Functional annotation]
    REF --> SnpEff
    
    SnpEff --> DbSNP[annotate_dbsnp<br/>rsID annotation]
    REF --> DbSNP
    
    DbSNP --> ClinVar[annotate_clinvar<br/>Clinical significance]
    REF --> ClinVar
    
    ClinVar --> COSMIC[annotate_cosmic<br/>Cancer mutations]
    REF --> COSMIC
    
    COSMIC --> FillTags[fill_tags<br/>Recalculate AC,AN,AF]
    
    FillTags --> Concat[concat_chunks_vcf<br/>Merge chunks per chr]
    
    Concat --> Singleton[singleton_dataset_chrm<br/>Identify singletons]
    Concat --> ExtractPGx[get_gene_vcf<br/>Extract PGx variants]
    
    ExtractPGx --> ConcatPGx[concat_vcf<br/>Combine chromosomes]
    ConcatPGx --> SitesOnly[sites_only<br/>Remove genotypes]
    ConcatPGx --> GetFreq[get_freq_3_2<br/>Calculate frequencies]
    
    SitesOnly --> MergeDS[merge_groups<br/>Merge all datasets]
    MergeDS --> FreqAll[get_freq_from_vcf_sites<br/>Calculate frequencies]
    FreqAll --> AnnotCV[add_clinvar_to_freq<br/>Add PharmGKB, GWAS]
    
    subgraph "Outputs"
        OUT1[Annotated VCFs]
        OUT2[PGx Variants]
        OUT3[Singleton Sites]
        OUT4[Merged Dataset]
    end
    
    Concat --> OUT1
    ConcatPGx --> OUT2
    Singleton --> OUT3
    AnnotCV --> OUT4
    
    style QC fill:#ffcccc
    style SnpEff fill:#ccffcc
    style DbSNP fill:#ccffcc
    style ClinVar fill:#ccffcc
    style COSMIC fill:#ccffcc
```

## PHASE 2: Population Grouping

```mermaid
graph TB
    subgraph "data_proc_pop"
        AnnVCF[Annotated VCFs<br/>from Phase 1] --> SplitPop[split_dataset_vcf_pop<br/>Split by population]
        SplitPop --> ExtractPGx[get_pgx_pop<br/>Extract PGx per pop]
        ExtractPGx --> ConcatPop[concat_pop<br/>Concat chromosomes]
    end
    
    subgraph "Grouping Workflows"
        ConcatPop --> GroupChrm[group_pops_by_chrm<br/>Group by chromosome]
        ConcatPop --> GroupAll[group_pops<br/>Concat all chr]
        ConcatPop --> Group2Pop[group_2pops<br/>Pairwise combinations]
        
        ConcatPop --> GroupPGx[group_pops_pgx<br/>Group PGx variants]
        ConcatPop --> Group2PGx[group_2pops_pgx<br/>Pairwise PGx]
    end
    
    subgraph "Outputs"
        OUT1[Population VCFs<br/>per chromosome]
        OUT2[Population VCFs<br/>all chromosomes]
        OUT3[Pairwise Population<br/>Combinations]
    end
    
    GroupChrm --> OUT1
    GroupAll --> OUT2
    Group2Pop --> OUT3
    Group2PGx --> OUT3
    
    style SplitPop fill:#fff4cc
    style ExtractPGx fill:#fff4cc
    style ConcatPop fill:#fff4cc
```

## PHASE 3: Statistical Analyses

```mermaid
graph TB
    subgraph "Fisher's Exact Test"
        GroupPGx[Grouped PGx<br/>Populations] --> FisherGroup[fisher_group_pgx<br/>Test group differences]
        FisherGroup --> VCFtoPlink1[vcf_to_plink<br/>Convert format]
        VCFtoPlink1 --> GenPheno1[generate_pheno_fam_sample<br/>Create phenotype files]
        GenPheno1 --> FisherTest1[fisher_test_plink<br/>Run Fisher test]
        
        Pairs[Pairwise<br/>Populations] --> Fisher2Pop[fisher_2pops_pgx<br/>Pairwise Fisher]
        Fisher2Pop --> VCFtoPlink2[vcf_to_plink]
        VCFtoPlink2 --> GenPheno2[generate_pheno_fam_sample1]
        GenPheno2 --> FisherTest2[fisher_test_plink_pops]
        FisherTest2 --> FilterFisher[filter_fisher_test<br/>Filter results]
        FilterFisher --> CombineCSV[combine_csv]
    end
    
    subgraph "FST Analysis"
        Pairs --> FSTAnalysis[fst_analysis<br/>Calculate FST]
        FSTAnalysis --> GetWeir[get_fst_weir_estimates<br/>Extract estimates]
        GetWeir --> CombineFST[combine_fst_weir_estimates<br/>Combine results]
        GetWeir --> GetWeirCut[get_fst_weir_estimates_cutoff<br/>Filter by threshold]
        GetWeirCut --> CombineFSTCut[combine_fst_weir_estimates_cutoff]
        
        GetWeir --> MeanFST1[combine_fst_mean_weir_estimates_1<br/>Per population mean]
        GetWeir --> MeanFST[combine_fst_mean_weir_estimates<br/>Overall mean]
        
        MeanFST --> GenMatrix[generate_fst_matrix<br/>Create FST matrix]
    end
    
    subgraph "Outputs"
        OUT1[Fisher Test<br/>P-values]
        OUT2[FST Estimates<br/>Per pair]
        OUT3[FST Matrix<br/>All populations]
    end
    
    CombineCSV --> OUT1
    CombineFST --> OUT2
    GenMatrix --> OUT3
    
    style FisherTest1 fill:#e1f5ff
    style FisherTest2 fill:#e1f5ff
    style FSTAnalysis fill:#ffe1f5
```

## PHASE 4: HDV Detection (Highly Differentiated Variants)

```mermaid
graph TB
    PopPGx[Population<br/>PGx VCFs] --> GetFreq[get_freq_3_3<br/>Calculate frequencies]
    FisherP[Fisher Test<br/>P-values] --> Combine1{Combine}
    FSTP[FST Estimates<br/>Per population] --> Combine1
    
    GetFreq --> SplitID[split_id_freq_2<br/>Split IDs]
    SplitID --> Combine1
    
    Combine1 --> HDV[hdv_dataset<br/>Identify HDVs<br/>Multiple thresholds:<br/>fold=2, ac=3<br/>group_threshold=1-5]
    
    subgraph "combine_hdvs_pops: Aggregate HDVs"
        HDV --> CombAll[combine_hdvs<br/>Combine all HDVs]
        HDV --> CombBase[combine_hdvs_base<br/>Combine base HDVs]
        
        CombAll --> AggAll[combine_hdvs_all<br/>Aggregate across pops]
        CombBase --> AggBase[combine_hdvs_base_all<br/>Aggregate base]
    end
    
    subgraph "Outputs"
        OUT1[HDV Variants<br/>Per population]
        OUT2[HDV Variants<br/>All populations]
        OUT3[HDV Variants<br/>Base + Aggregated]
    end
    
    CombAll --> OUT1
    AggAll --> OUT2
    AggBase --> OUT3
    
    style HDV fill:#ccffcc
    style CombAll fill:#ffccff
    style AggAll fill:#ffccff
```

## PHASE 5: Consequence Analysis (Optional)

```mermaid
graph TB
    PGxData[PGx Dataset<br/>from Phase 1] --> Filter{Filter<br/>AIBST only}
    SOTerms[SO Terms File] --> Filter
    
    Filter --> CSQ[csq<br/>Predict consequences<br/>COMMENTED OUT]
    CSQ -.-> PlotCSQ[plot_csq<br/>Visualize consequences<br/>COMMENTED OUT]
    
    subgraph "Optional Workflows - Not Active"
        PCA[pca<br/>Population structure]
        CountAll[count_all_dataset<br/>Count all variants]
        CountPGx[count_pgx_dataset<br/>Count PGx variants]
    end
    
    style CSQ fill:#cccccc,stroke-dasharray: 5 5
    style PlotCSQ fill:#cccccc,stroke-dasharray: 5 5
    style PCA fill:#cccccc,stroke-dasharray: 5 5
    style CountAll fill:#cccccc,stroke-dasharray: 5 5
    style CountPGx fill:#cccccc,stroke-dasharray: 5 5
```

## Data Flow Summary

```mermaid
graph LR
    A[Raw VCFs] --> B[QC + Annotation]
    B --> C[Chunked Processing]
    C --> D[Population Split]
    D --> E1[Fisher Tests]
    D --> E2[FST Analysis]
    E1 --> F[HDV Detection]
    E2 --> F
    F --> G[Results]
    
    B -.-> H[Optional: CSQ]
    B -.-> I[Optional: PCA]
    B -.-> J[Optional: Counts]
    
    style A fill:#e1f5ff
    style B fill:#ccffcc
    style C fill:#fff4cc
    style D fill:#ffe1cc
    style E1 fill:#e1f5ff
    style E2 fill:#ffe1f5
    style F fill:#ccffcc
    style G fill:#ccffff
    style H fill:#cccccc,stroke-dasharray: 5 5
    style I fill:#cccccc,stroke-dasharray: 5 5
    style J fill:#cccccc,stroke-dasharray: 5 5
```

## Key Process Characteristics

### Resource-Intensive Processes (BigMem/ExtraBig)

- `annotate_snpeff` - Functional annotation
- `get_map` - Chromosome extraction
- `fst_analysis` - Population differentiation
- `merge_pop_groups` - Large VCF merging

### Parallel Processing

- Chunks: 25MB regions processed in parallel
- Chromosomes: All requested chromosomes processed simultaneously
- Populations: Per-population analyses run in parallel
- Pairwise: All population pairs processed concurrently

### Caching Strategy

All processes support Nextflow caching with `-resume`:

- Previously successful processes reuse cached results
- Only failed or new processes re-execute
- Dramatically reduces re-run time (hours → seconds)

## Legend

- 🟦 Blue: Input/Data processing
- 🟩 Green: Annotation/Analysis
- 🟨 Yellow: Transformation/Grouping
- 🟪 Purple: Statistical tests
- ⬜ Gray (dashed): Optional/Commented out
