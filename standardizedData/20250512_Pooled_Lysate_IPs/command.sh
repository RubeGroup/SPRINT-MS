#!/bin/bash
python3 build.py ../../rawData/20250512_Pooled_Lysate_IPs/20250512_Pooled_Lysate_IPs_DDA_MQ_noMatchBwRuns_LFQ_iBAQ_proteinGroups_filt.txt --json config.json --key intensity     --output 20250512_Pooled_Lysate_IPs.intensity.tsv
python3 build.py ../../rawData/20250512_Pooled_Lysate_IPs/20250512_Pooled_Lysate_IPs_DDA_MQ_noMatchBwRuns_LFQ_iBAQ_proteinGroups_filt.txt --json config.json --key lfq           --output 20250512_Pooled_Lysate_IPs.lfq.tsv
python3 build.py ../../rawData/20250512_Pooled_Lysate_IPs/20250507_Pooled_Lysates_IP_DIA_automaticCrossRunNorm_noImp_Report_iBAQ.tsv      --json config.json --key DIA-intensity --output 20250512_Pooled_Lysate_IPs.DIA-intensity.tsv
#python3 build.py ../../rawData/20250512_Pooled_Lysate_IPs/20250507_Pooled_Lysates_IP_DIA_automaticCrossRunNorm_noImp_Report_iBAQ.tsv      --json config.json --key DIA-intensity --output /dev/stdout | sed 's/\t\t/\t0.0\t/g' | sed 's/\t\t/\t0.0\t/g' > 20250512_Pooled_Lysate_IPs.DIA-intensity.tsv

