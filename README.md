# 🧬 Genome Analysis and Annotation Tool

This repository contains a **Python-based tool** for fetching, analyzing, annotating, and visualizing genomic sequences using **Biopython** and **Matplotlib**.  

The script allows users to:  
- Fetch sequences from **NCBI** using accession numbers  
- Perform sequence operations (transcription, translation, reverse transcription)  
- Analyze **GC content** and sequence statistics  
- Align sequences using **PairwiseAligner**  
- Annotate genomes with features  
- Visualize outputs with plots and summaries  

---

## ✨ Features

- 📥 **Fetch Genomic Sequences** directly from NCBI (GenBank format with annotations).  
- 📊 **Sequence Analysis**: Calculate sequence length and GC content.  
- 🔬 **Molecular Operations**: Transcription, translation, and reverse transcription.  
- 🔗 **Sequence Alignment**: Global pairwise alignment with Biopython’s `PairwiseAligner`.  
- 🧬 **Genome Annotation**: Add example features (gene, exon, CDS) if none exist.  
- 📂 **Save Annotated Genome** as a GenBank file.  
- 📈 **Visualization**:  
  - Sequence length (bar chart)  
  - GC content (pie chart)  
  - Genome browser–style feature plots  
- 📜 **Text-based Outputs** for lightweight summaries.  

---

## 🛠 Requirements

Install dependencies before running:

```bash
pip install biopython matplotlib numpy
