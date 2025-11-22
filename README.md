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

This method has been validated on *Staphylococcus aureus* ST strains and extended toward *Mycobacterium tuberculosis* drug-resistance gene panels. :contentReference[oaicite:1]{index=1}

---

## 🚀 Features

### **1. Golden-Gate–Aided Within-Read Barcoding**
- Embeds unique barcodes *inside the read* rather than at adapters.
- Enables simultaneous assembly of 7–10 target amplicons into a single long fragment.  
- Verified feasibility via electrophoresis and long-read assembly. :contentReference[oaicite:2]{index=2}  
- Supports ≥108 samples in a single run.

### **2. Automated Liquid Handling**
- Supports OT-2 layout-to-layout coordination for reproducible multiplex PCR, ligation, purification, and pPCR workflows.
- Human-machine shared CSV layout interface improves transparency, reduces errors, and accelerates setup.  
  (See project slide deck.) :contentReference[oaicite:3]{index=3}

### **3. Barcode Compatibility & Optimization Tools**
- R scripts for primer overhang design, barcode compatibility scoring, and ligation efficiency prediction.
- Code structured for:
  - evaluating BsaI-HFv2 overhang combinations,  
  - Golden-Gate ligation models,  
  - multiplex PCR primer ratio optimization.

### **4. Bioinformatics Pipeline**
Includes:
- read trimming & barcode detection  
- barcode demultiplexing for within-read barcodes  
- error tolerance handling for Nanopore signals  



---

## 📁 Repository Structure

