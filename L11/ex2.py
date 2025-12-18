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

def find_detailed_connections(covid_seq, flu_seq, window_size=150, step=150, threshold=55):
    """
    Breaks COVID into windows, runs real local alignment against Flu,
    and stores the full aligned sequences for the best hits.
    """
    hits = []
    total_windows = len(covid_seq) // step
    print(f"Scanning {total_windows} windows (Size: {window_size}, Step: {step})...")
    print("Performing full local alignments (this takes time)...")
    
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

                    hits.append({
                        'c_start': start_c,
                        'c_end': end_c,
                        'f_start': aln.start, 
                        'f_end': aln.end,     
                        'score': score,
                        'aln_seq_c': aln.seqA[aln.start:aln.end],
                        'aln_seq_f': aln.seqB[aln.start:aln.end]
                    })
        except Exception as e:
            print(f"\nWarning in alignment window {i}: {e}")
            continue
            
    print("\nAlignment complete.")
    hits.sort(key=lambda x: x['score'], reverse=True)
    return hits[:8] 


def show_detail_plot(match_data, match_idx):
    """
    Pops up a new figure showing the base-by-base alignment text.
    """
    seq1 = match_data['aln_seq_c']
    seq2 = match_data['aln_seq_f']
    length = len(seq1)
    
    fig = plt.figure(figsize=(max(10, length/10), 4)) 
    ax = plt.gca()
    
    font_args = {'family': 'monospace', 'fontsize': 12, 'ha': 'center', 'va': 'center'}
    match_color = '#d4edda'
    mismatch_color = '#f8d7da'
    
    print(f"\n--- Detail View for Link {match_idx+1} ---")
    print(f"COVID: {seq1}")
    print(f"FLU  : {seq2}")

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

        rect = Rectangle((i - 0.5, 0), 1, 3, facecolor=bg_color, edgecolor='none', zorder=1)
        ax.add_patch(rect)

        plt.text(i, 2.2, base1, fontweight='bold', **font_args)
        plt.text(i, 1.5, connection_char, color='gray', **font_args)
        plt.text(i, 0.8, base2, fontweight='bold', **font_args)

    ax.set_yticks([0.8, 2.2])
    ax.set_yticklabels([f'Flu (Start: {match_data["f_start"]})', 
                        f'CoV (Start: {match_data["c_start"]})'])
    ax.set_xlim(-0.5, length - 0.5)
    ax.set_ylim(0, 3)
    ax.set_xticks(np.arange(0, length, 10))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.title(f"Detailed Alignment View: Link #{match_idx+1} (Score: {match_data['score']:.1f})")
    plt.tight_layout()
    plt.show()


def plot_interactive_connections(matches, len_covid, len_flu):
    fig, ax = plt.subplots(figsize=(14, 8))
    
    plt.hlines(1, 0, len_covid, color='#1f77b4', linewidth=6, label='SARS-CoV-2 (Query)')
    plt.hlines(0, 0, len_flu, color='#ff7f0e', linewidth=12, label='Influenza A (Target)')
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(matches)))
    pickable_artists = [] 

    print(f"Plotting {len(matches)} connections. CLICK on lines/boxes to see details.")

    for i, match in enumerate(matches):
        c_pos = (match['c_start'] + match['c_end'])/2
        f_pos = (match['f_start'] + match['f_end'])/2
        color = colors[i]
        
        line, = plt.plot([c_pos, f_pos], [1, 0], color=color, linestyle='--', 
                         alpha=0.7, linewidth=2, picker=5)
    
        label_text = (f"LINK {i+1}\nCLICK TO ZOOM\n"
                      f"C:{match['c_start']} \u2194 F:{match['f_start']}")
        
        text_box = plt.text((c_pos+f_pos)/2, 0.5, label_text, fontsize=8, 
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.9),
                 ha='center', va='center', color='black', picker=True)
                 
        pickable_artists.append((line, i))
        pickable_artists.append((text_box, i))

        plt.plot([match['c_start'], match['c_end']], [1, 1], color='red', linewidth=8)
        plt.plot([match['f_start'], match['f_end']], [0, 0], color='red', linewidth=14)

    def on_pick(event):
        artist = event.artist
        for pickable, idx in pickable_artists:
            if artist == pickable:
                show_detail_plot(matches[idx], idx)
                break

    fig.canvas.mpl_connect('pick_event', on_pick)

    plt.yticks([0, 1], ['Influenza A (PB2)', 'SARS-CoV-2'], fontsize=12, fontweight='bold')
    plt.xlabel('Base Pairs (bp)')
    plt.title(f'INTERACTIVE MAP: Top {len(matches)} Local Alignments\n(Click connecting lines or boxes to view sequence details)', fontsize=14, color='blue')
    plt.xlim(-500, max(len_covid, len_flu)+500)
    plt.ylim(-0.3, 1.3)
    plt.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    plt.show()

real_matches = find_detailed_connections(seq_covid, seq_flu, window_size=200, step=200, threshold=55)

if real_matches:
    plot_interactive_connections(real_matches, len(seq_covid), len(seq_flu))
else:
    print("\nNo significant biological alignments found above threshold.")
    print("(This is expected, as SARS-CoV-2 and Influenza are not genetically related).")
    print("Try lowering the threshold in 'find_detailed_connections' to see lower quality matches.")