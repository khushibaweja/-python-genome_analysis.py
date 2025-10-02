# 🧬 Genome Analysis and Annotation Tool  

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)  
![Biopython](https://img.shields.io/badge/Biopython-1.81%2B-green)  
![Matplotlib](https://img.shields.io/badge/Matplotlib-3.7%2B-orange)  
![NumPy](https://img.shields.io/badge/NumPy-1.24%2B-lightblue)  
![License](https://img.shields.io/badge/License-MIT-yellow)  

A **Python-based bioinformatics tool** for fetching, analyzing, annotating, and visualizing genomic sequences using **Biopython**, **Matplotlib**, and **NumPy**.  

This tool is designed for students, researchers, and bioinformatics enthusiasts who want to:  
- Retrieve annotated sequences directly from **NCBI**  
- Perform **molecular operations** like transcription, translation, and reverse transcription  
- Analyze **GC content, sequence length, and features**  
- Run **pairwise sequence alignments**  
- Annotate and export genomes in **GenBank format**  
- Generate **biological data visualizations**  

---

## ✨ Features  

- 📥 **Fetch Sequences**: Retrieve annotated genomic sequences from NCBI (GenBank format).  
- 📊 **Sequence Analysis**: Compute sequence length and GC content.  
- 🔬 **Molecular Operations**: Transcription, translation, and reverse transcription.  
- 🔗 **Sequence Alignment**: Global pairwise alignment using Biopython’s `PairwiseAligner`.  
- 🧬 **Genome Annotation**: Add example features (gene, exon, CDS) if none exist.  
- 📂 **Save Output**: Write annotated genome to a GenBank file.  
- 📈 **Visualization**:  
  - Sequence length (bar chart)  
  - GC content (pie chart)  
  - Genome browser–style feature plots  
- 📜 **Text Summaries**: Lightweight console-based visualization of sequence statistics.  

---

## 🛠 Tech Stack  

- **Python 3.8+**  
- [Biopython](https://biopython.org/) (`Bio.Entrez`, `SeqIO`, `SeqFeature`, `PairwiseAligner`)  
- [Matplotlib](https://matplotlib.org/) (for plotting)  
- [NumPy](https://numpy.org/) (for numerical operations)  

---

## 📦 Installation  

Clone this repository and install the required dependencies:  

```bash
git clone https://github.com/khushibaweja/genome-analysis-tool.git
cd genome-analysis-tool
