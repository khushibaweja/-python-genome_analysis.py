from Bio import Entrez, SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt
import numpy as np

# Set your email for NCBI compliance
Entrez.email = "kh13042004@gmail.com"  # Replace with your email

def fetch_sequence_from_accession(acc):
    """Fetch a nucleotide sequence from NCBI by accession number in GenBank format to include annotations."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        print(f"Fetched sequence {record.id} from NCBI with annotations.")
        print(f"Number of features: {len(record.features) if hasattr(record, 'features') else 0}")
        return record
    except Exception as e:
        print(f"Error fetching accession {acc}: {e}")
        return None

def analyze_features(record):
    """Calculate length and GC content of a sequence."""
    seq_str = str(record.seq).upper()
    length = len(seq_str)
    gc_count = seq_str.count("G") + seq_str.count("C")
    gc_content = 100 * gc_count / length if length > 0 else 0
    print(f"\n=== Sequence Features ===")
    print(f"Sequence length: {length} bp")
    print(f"GC content: {gc_content:.2f}%")
    if hasattr(record, 'features') and record.features:
        print(f"Number of annotations: {len(record.features)}")
    return length, gc_content

def safe_translate(seq):
    """Translate sequence safely by padding to multiple of 3."""
    remainder = len(seq) % 3
    if remainder != 0:
        seq += "N" * (3 - remainder)
    return seq.translate(to_stop=False)

def sequence_operations(dna_seq):
    """Perform transcription, translation, and reverse transcription."""
    rna_seq = dna_seq.transcribe()
    print("\n=== RNA Sequence ===")
    print(rna_seq)

    protein = safe_translate(rna_seq)
    print("\n=== Protein Translation ===")
    print(protein)

    rev_dna = rna_seq.back_transcribe()
    print("\n=== Reverse-Transcribed DNA ===")
    print(rev_dna)

    return rna_seq, protein, rev_dna

def align_sequences(seq1, seq2):
    """Perform global pairwise alignment and print results."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignments = aligner.align(seq1, seq2)
    best = alignments[0]

    print(f"\n=== Pairwise Alignment ===")
    print(f"Alignment score: {best.score:.2f}\n")

    # Access aligned sequences (Biopython 1.81+)
    seq1_aligned, seq2_aligned = best.sequences
    print("Sequence 1 (aligned):")
    print(seq1_aligned)
    print("\nSequence 2 (aligned):")
    print(seq2_aligned)

    return best

def annotate_genome(record):
    """Enhance annotations if needed (e.g., add example if no features exist)."""
    if not hasattr(record, 'features') or not record.features:
        print("No existing features found. Adding example gene and exon annotations.")
        features = [
            SeqFeature.SeqFeature(
                SeqFeature.FeatureLocation(0, 100, strand=1),
                type="gene",
                qualifiers={"note": "Example gene", "gene": "EXAMPLE"}
            ),
            SeqFeature.SeqFeature(
                SeqFeature.FeatureLocation(10, 50, strand=1),
                type="exon",
                qualifiers={"note": "Example exon"}
            ),
            SeqFeature.SeqFeature(
                SeqFeature.FeatureLocation(20, 40, strand=1),
                type="CDS",
                qualifiers={"note": "Example CDS", "product": "Example protein"}
            )
        ]
        if not hasattr(record, "features"):
            record.features = []
        record.features.extend(features)
        print("Added example gene, exon, and CDS annotations to the genome.")
    else:
        print(f"Genome already has {len(record.features)} annotations. No additions made.")
    return record

def plot_gc_and_length(length, gc_content):
    """Visualize sequence length and GC content using matplotlib."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Bar plot for sequence length
    ax1.bar(['Sequence Length (bp)'], [length], color='skyblue')
    ax1.set_ylabel('Length (bp)')
    ax1.set_title('Sequence Length')
    ax1.tick_params(axis='x', rotation=0)

    # Pie chart for GC content
    ax2.pie([gc_content, 100 - gc_content], labels=['GC (%)', 'AT (%)'], 
            autopct='%1.1f%%', colors=['darkgreen', 'orange'], startangle=90)
    ax2.set_title('GC Content Distribution')

    plt.tight_layout()
    plt.show()

def plot_features(record):
    """Plot gene annotations and features as a simple genome browser view using matplotlib."""
    if not hasattr(record, 'features') or not record.features:
        print("No features to plot.")
        return

    fig, ax = plt.subplots(figsize=(12, max(4, len(record.features) * 0.8)))
    length = len(record.seq)
    ax.set_xlim(0, length)
    ax.set_ylim(0, len(record.features) + 1)
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Features')
    ax.set_title(f'Gene Annotations and Features in {record.id}')

    colors = {'gene': 'blue', 'exon': 'green', 'CDS': 'red', 'mRNA': 'purple', 'five_prime_UTR': 'yellow', 'three_prime_UTR': 'orange'}
    y_positions = np.arange(1, len(record.features) + 1)

    for i, feature in enumerate(record.features):
        start = feature.location.start.position if hasattr(feature.location, 'position') else feature.location.start
        end = feature.location.end.position if hasattr(feature.location, 'position') else feature.location.end
        feat_type = feature.type
        color = colors.get(feat_type, 'gray')
        
        # Draw rectangle for feature
        ax.add_patch(plt.Rectangle((start, i + 0.2), end - start, 0.6, 
                                   facecolor=color, alpha=0.7, edgecolor='black'))
        
        # Add label
        ax.text((start + end) / 2, i + 0.5, f"{feat_type}", ha='center', va='center', 
                fontsize=8, rotation=0, bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))

    # Draw the sequence backbone
    ax.plot([0, length], [0.5, 0.5], 'k-', linewidth=2, label='Sequence Backbone')
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

def visualize_outputs_text(length, gc_content):
    """Text-based visualization of sequence length and GC content."""
    print("\n=== Visual Summary (Text) ===")
    
    # Sequence length visualization
    print(f"\nSequence Length: {length} bp")
    scale_len = max(length // 50, 1)
    print("Length: " + "*" * (length // scale_len))
    
    # GC content visualization
    print(f"\nGC Content: {gc_content:.2f}%")
    scale_gc = max(int(gc_content // 2), 1)
    print("GC: " + "*" * scale_gc)

def main():
    accession = input("Enter the accession number: ").strip() or "NM_000249"
    record = fetch_sequence_from_accession(accession)
    if not record:
        print("Failed to fetch sequence. Exiting.")
        return

    # 1. Analyze features
    length, gc_content = analyze_features(record)

    # 2. Sequence operations (transcription, translation, reverse transcription)
    rna_seq, protein_seq, rev_dna = sequence_operations(record.seq)

    # 3. Align DNA sequence to itself (example alignment)
    alignment = align_sequences(record.seq, record.seq)

    # 4. Annotate genome (enhance if needed)
    annotated_record = annotate_genome(record)

    # 5. Add molecule_type to avoid GenBank write error if missing
    if "molecule_type" not in annotated_record.annotations:
        annotated_record.annotations["molecule_type"] = "DNA"

    # 6. Save annotated genome to GenBank file
    genbank_filename = f"annotated_{record.id}.gb"
    SeqIO.write(annotated_record, genbank_filename, "genbank")
    print(f"\nAnnotated genome saved to {genbank_filename}")

    # 7. Visualize outputs as text
    visualize_outputs_text(length, gc_content)

    # 8. Matplotlib visualizations
    print("\nGenerating matplotlib plots...")
    plot_gc_and_length(length, gc_content)
    plot_features(annotated_record)

# ✅ Corrected entry point
if __name__ == "__main__":
    main()
