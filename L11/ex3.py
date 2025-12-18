import numpy as np
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO, pairwise2
from matplotlib.patches import Rectangle
import os

Entrez.email = "vlad.paraschiv1511@stud.fils.upb.ro" 

def download_genome(accession_id, filename):
    if os.path.exists(filename):
        print(f"Loading {filename} locally...")
        record = SeqIO.read(filename, "fasta")
    else:
        print(f"Downloading {accession_id}...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            with open(filename, "w") as output_handle:
                SeqIO.write(record, output_handle, "fasta")
        except Exception as e:
            print(f"Error: {e}")
            return ""
    return str(record.seq).upper()

print("--- Fetching Genomes ---")
seq_covid = download_genome("NC_045512", "covid19.fasta")
seq_flu = download_genome("NC_002019", "flu_segment1.fasta")

if not seq_covid or not seq_flu:
    print("Failed to load genomes. Exiting.")
    exit()

class SimilarityScorer:
    """Three quantitative similarity scoring equations for aligned sequences."""
    
    @staticmethod
    def score_identity(align1: str, align2: str) -> tuple[float, str]:
        """1. IDENTITY SCORE: % identical positions (excluding gaps)"""
        matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
        total_positions = sum(1 for a, b in zip(align1, align2) if a != '-' and b != '-')
        identity = (matches / total_positions * 100) if total_positions > 0 else 0
        return identity, f"{identity:.1f}% ID ({matches}/{total_positions})"
    
    @staticmethod
    def score_similarity_matrix(align1: str, align2: str) -> tuple[float, str]:
        """2. SIMILARITY MATRIX: BLOSUM62-inspired nucleotide scoring"""
        nuc_sim = {
            ('A','A'): 5,  ('T','T'): 5,  ('C','C'): 5,  ('G','G'): 5,
            ('A','G'): 2,  ('G','A'): 2,  ('C','T'): 2,  ('T','C'): 2,  # Transitions
            ('A','T'): -3, ('T','A'): -3, ('A','C'): -3, ('C','A'): -3,
            ('G','T'): -3, ('T','G'): -3, ('G','C'): -3, ('C','G'): -3
        }
        
        total_score, count = 0, 0
        for a, b in zip(align1, align2):
            if a == '-' or b == '-': continue
            score = nuc_sim.get((a, b), -4)
            total_score += score
            count += 1
        
        avg_similarity = total_score / count if count > 0 else 0
        return avg_similarity, f"{avg_similarity:.1f} avg sim ({count} pos)"
    
    @staticmethod
    def score_conservation(align1: str, align2: str) -> tuple[float, str]:
        """3. CONSERVATION SCORE: Identity + gap balance penalty"""
        matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
        aligned_len = len(align1)
        gaps1, gaps2 = align1.count('-'), align2.count('-')
        gap_penalty = abs(gaps1 - gaps2) * 0.5
        
        conservation = (matches / max(1, aligned_len - gap_penalty)) * 100
        return conservation, f"{conservation:.1f}% cons (gap_pen: {gap_penalty:.0f})"

def find_detailed_connections(covid_seq, flu_seq, window_size=150, step=150, threshold=55):
    """
    Enhanced version with THREE similarity scoring equations applied to each alignment.
    """
    hits = []
    scorer = SimilarityScorer()
    total_windows = len(covid_seq) // step
    print(f"Scanning {total_windows} windows (Size: {window_size}, Step: {step})...")
    print("Performing full local alignments + 3 similarity scores...")
    
    for i, start_c in enumerate(range(0, len(covid_seq) - window_size, step)):
        if i % 5 == 0: print(f"Processing window {i}/{total_windows}...", end='\r')
        
        chunk_c = covid_seq[start_c : start_c + window_size]

        try:
            alignments = pairwise2.align.localms(chunk_c, flu_seq, 3, -2, -3, -1, one_alignment_only=True)
            
            if alignments:
                aln = alignments[0]
                score = aln.score
                
                if score > threshold:
                    ungapped_len_c = len(aln.seqA.replace('-', ''))
                    end_c = start_c + ungapped_len_c

                    # === APPLY THREE SCORING EQUATIONS ===
                    id_score, id_desc = scorer.score_identity(aln.seqA, aln.seqB)
                    sim_score, sim_desc = scorer.score_similarity_matrix(aln.seqA, aln.seqB)
                    cons_score, cons_desc = scorer.score_conservation(aln.seqA, aln.seqB)
                    
                    hits.append({
                        'c_start': start_c,
                        'c_end': end_c,
                        'f_start': aln.start, 
                        'f_end': aln.end,     
                        'raw_score': score,
                        'id_score': id_score,
                        'sim_score': sim_score,
                        'cons_score': cons_score,
                        'id_desc': id_desc,
                        'sim_desc': sim_desc,
                        'cons_desc': cons_desc,
                        'aln_seq_c': aln.seqA,
                        'aln_seq_f': aln.seqB
                    })
        except Exception as e:
            print(f"\nWarning in window {i}: {e}")
            continue
            
    print("\nAlignment + scoring complete.")
    # Sort by identity score (most biologically meaningful)
    hits.sort(key=lambda x: x['id_score'], reverse=True)
    return hits[:8] 

def show_detail_plot(match_data, match_idx):
    """
    Enhanced detail view showing all THREE scores.
    """
    seq1 = match_data['aln_seq_c']
    seq2 = match_data['aln_seq_f']
    length = len(seq1)
    
    fig = plt.figure(figsize=(max(12, length/8), 5)) 
    ax = plt.gca()
    
    font_args = {'family': 'monospace', 'fontsize': 11, 'ha': 'center', 'va': 'center'}
    match_color = '#d4edda'
    mismatch_color = '#f8d7da'
    
    print(f"\n--- Detail View for Link {match_idx+1} ---")
    print(f"COVID: {seq1}")
    print(f"FLU:   {seq2}")
    print(f"Scores: RAW={match_data['raw_score']:.1f} | {match_data['id_desc']} | {match_data['sim_desc']} | {match_data['cons_desc']}")

    # Three scoring equations summary
    print("\n=== SCORING EQUATIONS ===")
    print("1. IDENTITY:     Matches / Total_positions × 100")
    print("2. SIMILARITY:   Σ(nuc_matrix_scores) / N")
    print("3. CONSERVATION: Matches / (Len - Gap_penalty) × 100")

    for i in range(length):
        base1 = seq1[i]
        base2 = seq2[i]
       
        if base1 == base2 and base1 not in ['-', 'N']:
            bg_color = match_color
            connection_char = "|"
        elif base1 == '-' or base2 == '-':
            bg_color = 'white'
            connection_char = " "
        else:
            bg_color = mismatch_color
            connection_char = "."

        rect = Rectangle((i - 0.5, 0), 1, 3.5, facecolor=bg_color, edgecolor='none', zorder=1)
        ax.add_patch(rect)

        plt.text(i, 2.6, base1, fontweight='bold', **font_args)
        plt.text(i, 1.8, connection_char, color='gray', **font_args)
        plt.text(i, 1.0, base2, fontweight='bold', **font_args)

    ax.set_yticks([1.0, 2.6])
    ax.set_yticklabels([f'Flu PB2 (Start: {match_data["f_start"]})', 
                        f'CoV-2 (Start: {match_data["c_start"]})'])
    ax.set_xlim(-0.5, length - 0.5)
    ax.set_ylim(0, 3.5)
    ax.set_xticks(np.arange(0, length, 10))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add scoring summary box
    scores_text = (f"RAW: {match_data['raw_score']:.1f}\n"
                   f"ID: {match_data['id_desc']}\n"
                   f"SIM: {match_data['sim_desc']}\n"
                   f"CONS: {match_data['cons_desc']}")
    plt.text(0.02, 0.98, scores_text, transform=ax.transAxes, 
             bbox=dict(boxstyle="round,pad=0.5", fc="lightblue", ec="blue", alpha=0.8),
             va='top', fontsize=10, fontweight='bold')
    
    plt.title(f"Detailed Alignment + 3 Similarity Scores: Link #{match_idx+1}", fontsize=14)
    plt.tight_layout()
    plt.show()

def plot_interactive_connections(matches, len_covid, len_flu):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), 
                                   gridspec_kw={'height_ratios': [3, 1]})
    
    # Main alignment plot
    ax1.hlines(1, 0, len_covid, color='#1f77b4', linewidth=6, label='SARS-CoV-2')
    ax1.hlines(0, 0, len_flu, color='#ff7f0e', linewidth=12, label='Influenza A PB2')
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(matches)))
    pickable_artists = [] 

    print(f"Plotting {len(matches)} scored connections. CLICK lines/boxes for details.")

    for i, match in enumerate(matches):
        c_pos = (match['c_start'] + match['c_end'])/2
        f_pos = (match['f_start'] + match['f_end'])/2
        color = colors[i]
        
        line, = ax1.plot([c_pos, f_pos], [1, 0], color=color, linestyle='--', 
                         alpha=0.7, linewidth=3, picker=5)
        
        # Enhanced label with all 3 scores
        label_text = (f"#{i+1} ID:{match['id_score']:.0f}%\n"
                      f"C:{match['c_start']}-{match['c_end']} ↔ F:{match['f_start']}")
        
        text_box = ax1.text((c_pos+f_pos)/2, 0.5, label_text, fontsize=9, 
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.9),
                 ha='center', va='center', picker=True)
                 
        pickable_artists.append((line, i))
        pickable_artists.append((text_box, i))

        ax1.plot([match['c_start'], match['c_end']], [1, 1], color=color, linewidth=10, alpha=0.6)
        ax1.plot([match['f_start'], match['f_end']], [0, 0], color=color, linewidth=16, alpha=0.6)

    # Scores comparison bar chart
    match_ids = [m['id_score'] for m in matches[:8]]
    match_sims = [m['sim_score'] for m in matches[:8]]
    match_cons = [m['cons_score'] for m in matches[:8]]
    
    x = np.arange(len(matches[:8]))
    width = 0.25
    
    ax2.bar(x - width, match_ids[:8], width, label='Identity %', alpha=0.8, color='green')
    ax2.bar(x, match_sims[:8], width, label='Similarity Avg', alpha=0.8, color='orange')
    ax2.bar(x + width, match_cons[:8], width, label='Conservation %', alpha=0.8, color='purple')
    
    ax2.set_xlabel('Top Matches')
    ax2.set_ylabel('Score Value')
    ax2.set_title('THREE Similarity Scoring Equations Comparison')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    # Main plot formatting
    ax1.set_yticks([0, 1])
    ax1.set_yticklabels(['Influenza A PB2 (2341 bp)', 'SARS-CoV-2 (29,903 bp)'], fontweight='bold')
    ax1.set_xlabel('Position (bp)')
    ax1.set_title(f'Top {len(matches)} Local Alignments with 3 Similarity Scores\nClick connections for detailed sequence + scoring analysis', fontsize=14)
    ax1.set_xlim(-500, max(len_covid, len_flu)+500)
    ax1.set_ylim(-0.3, 1.3)
    ax1.grid(axis='x', alpha=0.3)
    ax1.legend()

    def on_pick(event):
        artist = event.artist
        for pickable, idx in pickable_artists:
            if artist == pickable:
                show_detail_plot(matches[idx], idx)
                break

    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.tight_layout()
    plt.show()

# Execute with enhanced scoring
real_matches = find_detailed_connections(seq_covid, seq_flu, window_size=200, step=200, threshold=45)

if real_matches:
    print("\n=== TOP MATCHES WITH THREE SCORING EQUATIONS ===")
    for i, match in enumerate(real_matches[:5]):
        print(f"{i+1}. ID:{match['id_score']:.1f}% | SIM:{match['sim_score']:.1f} | CONS:{match['cons_score']:.1f}% | RAW:{match['raw_score']:.1f}")
    
    plot_interactive_connections(real_matches, len(seq_covid), len(seq_flu))
else:
    print("No significant alignments found. Lower threshold in find_detailed_connections().")
