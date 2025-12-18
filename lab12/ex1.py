import math
import matplotlib.pyplot as plt
plt.switch_backend('TkAgg') 
from Bio import Entrez, SeqIO

Entrez.email = "bugamihai157@gmail.com"

motifs = [
    "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
    "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
]
motif_length = 9
null_model = 0.25

influenza_accessions = [
    "NC_026433", "NC_026434", "NC_026435", "NC_026436", "NC_026437",
    "NC_026438", "NC_026431", "NC_026432", "NC_002017", "NC_002018"
]

def calculate_log_matrix():
    counts = {'A': [0]*motif_length, 'C': [0]*motif_length, 'G': [0]*motif_length, 'T': [0]*motif_length}
    for seq in motifs:
        for i, nuc in enumerate(seq):
            counts[nuc][i] += 1
            
    total = len(motifs)
    logs = {'A': [], 'C': [], 'G': [], 'T': []}
    
    for nuc in ['A', 'C', 'G', 'T']:
        for i in range(motif_length):
            p = counts[nuc][i] / total
            if p > 0:
                score = math.log(p / null_model)
            else:
                score = -15.0
            logs[nuc].append(score)
    return logs

def scan_genome(sequence, log_matrix):
    scores = []
    positions = []
    seq_upper = str(sequence).upper()
    
    for i in range(len(seq_upper) - motif_length + 1):
        window = seq_upper[i : i + motif_length]
        current_score = 0
        valid_window = True
        
        for pos, nuc in enumerate(window):
            if nuc in log_matrix:
                current_score += log_matrix[nuc][pos]
            else:
                valid_window = False
                break
        
        positions.append(i)
        if valid_window:
            scores.append(current_score)
        else:
            scores.append(-15.0)
            
    return positions, scores

def run_analysis():
    log_matrix = calculate_log_matrix()
    
    fig, axes = plt.subplots(5, 2, figsize=(14, 12))
    axes = axes.flatten()
    fig.patch.set_facecolor('white')
    
    print(f"Downloading and scanning {len(influenza_accessions)} genomes...")
    
    for idx, accession in enumerate(influenza_accessions):
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            positions, scores = scan_genome(record.seq, log_matrix)
            
            ax = axes[idx]
            ax.set_facecolor('white')
            
            baseline = -15.0
            ax.vlines(positions, ymin=baseline, ymax=scores, color='#1f77b4', linewidth=1.2, zorder=1)
            ax.axhline(0, color='#d62728', linestyle='--', linewidth=1.0, zorder=2)
            
            high_signals_x = []
            high_signals_y = []
            max_score = -999
            max_pos = -1
            
            for p, s in zip(positions, scores):
                if s > 0:
                    high_signals_x.append(p)
                    high_signals_y.append(s)
                    if s > max_score:
                        max_score = s
                        max_pos = p
            
            ax.scatter(high_signals_x, high_signals_y, color='green', s=25, zorder=3, edgecolors='black', linewidth=0.5)
            
            if max_score > 0:
                ax.annotate(f"Hit: {max_pos}\n{max_score:.1f}", 
                            xy=(max_pos, max_score), 
                            xytext=(max_pos + 50, max_score + 3),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                            fontsize=8, ha='left')

            ax.set_title(f"{accession} ({len(record.seq)} bp)", fontsize=10, fontweight='bold')
            ax.set_ylabel("Score", fontsize=8)
            ax.set_ylim(bottom=-16, top=10) 
            
            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.grid(True, which='major', axis='y', linestyle=':', linewidth=0.5, color='gray', alpha=0.3)
            
            if idx >= 8:
                ax.set_xlabel("Genome Position (bp)", fontsize=9)

        except Exception as e:
            print(f"Error processing {accession}: {e}")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.2)
    print("Done. Displaying plots...")
    plt.show()

if __name__ == "__main__":
    run_analysis()