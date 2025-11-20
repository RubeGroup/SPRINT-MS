#!/binb/bash
python3 build.py ../../rawData/20250404_PooledIP_15x30/20250404_PooledIP_15x30_DDA_MQ_noMatchBwRuns_LFQ_iBAQ_proteinGroups.txt --json config.json --key intensity --output 20250404_PooledIP_15x30.intensity.tsv
python3 build.py ../../rawData/20250404_PooledIP_15x30/20250404_PooledIP_15x30_DDA_MQ_noMatchBwRuns_LFQ_iBAQ_proteinGroups.txt --json config.json --key lfq --output 20250404_PooledIP_15x30.lfq.tsv
python3 build.py ../../rawData/20250404_PooledIP_15x30/20250624_Pooled_Ab_IP_15x30_DIA_noCrossRunNorm_noImp_Report.tsv         --json config.json --key DIA-intensity --output 20250404_PooledIP_15x30.DIA-intensity.tsv
