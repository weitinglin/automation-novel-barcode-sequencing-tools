# Novel Within-Read Barcoded Strategy for Multiplexing Targeted Nanopore Sequencing

This repository hosts the source code, analysis scripts, and automation templates used in the project **“Novel Within-Reads Barcoded Strategy Combined with Automation for Multiplexing Targeted Nanopore Sequencing.”**

Our method combines:
- a **Golden-Gate–assisted barcoded library construction workflow**,  
- an **automated liquid-handling pipeline**,  
- **targeted long-read sequencing (Oxford Nanopore)**, and  
- **bioinformatics tools** for barcode extraction, demultiplexing, and variant calling, MLST subtyping

The goal is to enable **low-cost, high-throughput targeted sequencing** for applications such as:
- *MLST genotyping*,  
- *high-multiplex amplicon sequencing*
- *pilot workflows for automated sequencing laboratories*.

This method has been validated on *Staphylococcus aureus* ST strains

---

## 🚀 Features

### **0. Ovewview:Golden-Gate–Aided Within-Read Barcoding**
- Embeds unique barcodes *inside the read* rather than at adapters.
- Enables simultaneous assembly of 7–10 target amplicons into a single long fragment.  
- Verified feasibility via electrophoresis and long-read assembly.
- Supports ≥108 samples in a single run.

### **1. Automated Liquid Handling**
- Jupyter notebook as opentron server usage,
    - Supports OT-2 layout-to-layout coordination for reproducible multiplex PCR, ligation, purification workflows.
    - Human-machine shared CSV layout interface improves transparency, reduces errors, and accelerates setup.  
  

### **2. Barcode Compatibility & Optimization Tools**
- R scripts for primer overhang design, barcode compatibility scoring, and ligation efficiency prediction.
- Code structured for:
  - evaluating BsaI-HFv2 overhang combinations,  
  - Golden-Gate ligation models,  
  - multiplex PCR primer ratio optimization.

### **3. Bioinformatics Pipeline**
- Shellscripts based, bioinformatic software involved Seqkit, bioawk, flye, bwa, samtools, IGV, medaka, minimap2
    - read trimming & barcode detection  
    - barcode demultiplexing for within-read barcodes  
    - assembly contig from long reads data
    - error tolerance experiments for Nanopore signals  
    - MLST subtyping and visualization result preparation

### **4. Data layer interface reference files
- step 1 : 16tube_position.csv, 20220816_16_endproductLayout_set.csv, barcode108_combination.csv, DNAtoPosition.csv
- step 2 : 16_endproductLayout_set.csv,16new_endproductLayout_set.csv
- step 3 : 96plate_endproduct.csv, 96plate_layout.csv, 16_endproductLayout_setv1.csv, 16new_endproductLayout.csv
- step 4 : 38_input_layout.csv, 38_output_layout.csv



