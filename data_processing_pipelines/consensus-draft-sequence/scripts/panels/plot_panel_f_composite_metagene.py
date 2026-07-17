#!/usr/bin/env python3
import gzip
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import gaussian_kde

def get_length(intervals):
    return sum(e - s + 1 for s, e in intervals)

def get_relative_pos(pos, intervals, strand):
    total_len = get_length(intervals)
    if total_len == 0:
        return None
    sorted_ivs = sorted(intervals, reverse=(strand == '-'))
    dist = 0
    for s, e in sorted_ivs:
        if s <= pos <= e:
            if strand == '+':
                dist += (pos - s + 1)
            else:
                dist += (e - pos + 1)
            return (dist - 1) / total_len
        else:
            dist += (e - s + 1)
    return None

def main():
    from pathlib import Path
    GITHUB_ROOT = Path(__file__).resolve().parents[2]
    
    LOCAL_GTF = GITHUB_ROOT / "inputs" / "gencode.v49.primary_assembly.annotation.gtf.gz"
    FALLBACK_GTF = Path.home() / "ref" / "gencode.v49.primary_assembly.annotation.gtf.gz"
    gtf_path = str(LOCAL_GTF if LOCAL_GTF.exists() else FALLBACK_GTF)
    
    if not Path(gtf_path).exists():
        raise FileNotFoundError(f"GENCODE GTF not found at {LOCAL_GTF} or {FALLBACK_GTF}.")
        
    tsv_path = str(GITHUB_ROOT / "outputs" / "tiered_lists" / "tiered_polyA.tsv")
    outdir = str(GITHUB_ROOT / "figures" / "panel_f_composite_metagene")
    os.makedirs(outdir, exist_ok=True)
    
    print("Parsing GTF (filtering for protein-coding genes/transcripts)...")
    BIN_SIZE = 100000
    transcripts = {}
    
    with gzip.open(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if 'gene_type "protein_coding"' not in line:
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ['exon', 'CDS', 'UTR']:
                continue
            
            seqname = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            
            attr = parts[8]
            tid = None
            gid = None
            for a in attr.split(';'):
                a = a.strip()
                if a.startswith('transcript_id'):
                    tid = a.split(' ')[1].strip('"')
                elif a.startswith('gene_id'):
                    gid = a.split(' ')[1].strip('"')
                    
            if not tid or not gid:
                continue
            
            if tid not in transcripts:
                transcripts[tid] = {
                    'chr': seqname, 
                    'strand': strand, 
                    'gene_id': gid,
                    'exon': [], 
                    'CDS': [], 
                    'UTR': [], 
                    'five_prime_utr': [], 
                    'three_prime_utr': []
                }
                
            transcripts[tid][feature].append((start, end))

    print("Deducing UTRs and selecting representative transcripts...")
    for tid, t in transcripts.items():
        if t['CDS'] and t['UTR']:
            cds_min = min(s for s, e in t['CDS'])
            cds_max = max(e for s, e in t['CDS'])
            for s, e in t['UTR']:
                if e < cds_min:
                    if t['strand'] == '+':
                        t['five_prime_utr'].append((s, e))
                    else:
                        t['three_prime_utr'].append((s, e))
                elif s > cds_max:
                    if t['strand'] == '+':
                        t['three_prime_utr'].append((s, e))
                    else:
                        t['five_prime_utr'].append((s, e))

    gene_to_best_tx = {}
    for tid, tx in transcripts.items():
        if not tx['CDS'] or not tx['five_prime_utr'] or not tx['three_prime_utr'] or not tx['exon']:
            continue
        gene_id = tx['gene_id']
        cds_len = sum(e - s + 1 for s, e in tx['CDS'])
        if gene_id not in gene_to_best_tx or cds_len > gene_to_best_tx[gene_id]['cds_len']:
            gene_to_best_tx[gene_id] = {
                'tid': tid,
                'cds_len': cds_len,
                'tx': tx
            }
            
    representative_transcripts = {item['tid']: item['tx'] for item in gene_to_best_tx.values()}
    print(f"Selected {len(representative_transcripts)} representative transcripts.")
    
    bins = defaultdict(list)
    for tid, tx in representative_transcripts.items():
        min_s = min(s for s, e in tx['exon'])
        max_e = max(e for s, e in tx['exon'])
        tx['start'] = min_s
        tx['end'] = max_e
        
        start_bin = min_s // BIN_SIZE
        end_bin = max_e // BIN_SIZE
        for b in range(start_bin, end_bin + 1):
            bins[(tx['chr'], b)].append(tid)
            
    print(f"Reading modification sites...")
    sites = pd.read_csv(tsv_path, sep='\t')
    mod_types = sorted(sites['name'].unique())
    
    # Store mapped positions for each modification
    mapped_positions = {}
    
    for mod in mod_types:
        mod_sites = sites[sites['name'] == mod]
        metagene_positions = []
        
        for _, row in mod_sites.iterrows():
            pos = row['start'] + 1
            chrom = row['chr']
            strand = row['strand']
            
            b = pos // BIN_SIZE
            candidate_tids = bins.get((chrom, b), [])
            
            overlapping_tids = []
            for tid in candidate_tids:
                tx = representative_transcripts[tid]
                if tx['strand'] != strand:
                    continue
                if tx['start'] <= pos <= tx['end']:
                    overlapping_tids.append(tid)
                    
            if not overlapping_tids:
                continue
                
            for tid in overlapping_tids:
                tx = representative_transcripts[tid]
                if any(s <= pos <= e for s, e in tx['five_prime_utr']):
                    rel = get_relative_pos(pos, tx['five_prime_utr'], strand)
                    if rel is not None:
                        metagene_positions.append(rel)
                        break
                elif any(s <= pos <= e for s, e in tx['CDS']):
                    rel = get_relative_pos(pos, tx['CDS'], strand)
                    if rel is not None:
                        metagene_positions.append(1.0 + rel)
                        break
                elif any(s <= pos <= e for s, e in tx['three_prime_utr']):
                    rel = get_relative_pos(pos, tx['three_prime_utr'], strand)
                    if rel is not None:
                        metagene_positions.append(2.0 + rel)
                        break
                        
        mapped_positions[mod] = np.array(metagene_positions)
        print(f"  {mod}: mapped {len(metagene_positions)} / {len(mod_sites)} sites")

    # Colors for each mod
    mod_colors = {
        'm6A': {'line': '#721817', 'fill': '#721817'},
        'm5C': {'line': '#001427', 'fill': '#001427'},
        'I':   {'line': '#2D6E1E', 'fill': '#2D6E1E'},
        'Y':   {'line': '#F0A202', 'fill': '#F0A202'},
    }

    # ----------------- PLOT 1: Overlay (Single Panel) -----------------
    # Sized for 1/4 of A4 page (roughly 4.13" x 5.85"). A size of 4.5" x 4.0" is perfect.
    plt.figure(figsize=(4.5, 3.8))
    plt.grid(axis='y', linestyle=':', alpha=0.5, zorder=0)
    
    for mod in ["m6A", "m5C", "Y", "I"]:
        pos = mapped_positions[mod]
        if len(pos) == 0:
            continue
        try:
            kde = gaussian_kde(pos)
            x_eval = np.linspace(0, 3, 300)
            y_eval = kde(x_eval)
        except:
            counts, edges = np.histogram(pos, bins=100, range=(0, 3), density=True)
            x_eval = (edges[:-1] + edges[1:]) / 2
            y_eval = counts
            
        colors = mod_colors.get(mod, {'line': '#455A64', 'fill': '#CFD8DC'})
        plt.plot(x_eval, y_eval, color=colors['line'], lw=2, label=mod, zorder=3)
        plt.fill_between(x_eval, 0, y_eval, color=colors['fill'], alpha=0.15, zorder=2)
        
    plt.axvline(1.0, color='#37474F', linestyle='--', lw=1.0, zorder=4)
    plt.axvline(2.0, color='#37474F', linestyle='--', lw=1.0, zorder=4)
    
    plt.xticks([0.5, 1.5, 2.5], ["5' UTR", "CDS", "3' UTR"], fontsize=10, fontweight='bold')
    plt.xlim(0, 3)
    plt.ylim(bottom=0)
    plt.xlabel("Transcript Region", fontsize=11, fontweight='bold', labelpad=6)
    plt.ylabel("Density", fontsize=11, fontweight='bold', labelpad=6)
    
    plt.legend(frameon=True, facecolor='white', edgecolor='none', fontsize=9, loc='upper left')
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#37474F')
    ax.spines['bottom'].set_color('#37474F')
    ax.tick_params(colors='#37474F', labelsize=9)
    
    plt.tight_layout()
    overlay_png = os.path.join(outdir, "panel_f_composite_metagene_overlay.png")
    overlay_pdf = os.path.join(outdir, "panel_f_composite_metagene_overlay.pdf")
    plt.savefig(overlay_png, dpi=300)
    plt.savefig(overlay_pdf)
    plt.close()
    print(f"Saved overlay plot to {overlay_png}")

    # ----------------- PLOT 2: 4x1 Stacked Column (1/3 A4 width) -----------------
    # 1/3 A4 width ≈ 70 mm ≈ 2.76 inches; tall stack for 4 modification types.
    fig, axes = plt.subplots(4, 1, figsize=(2.8, 7.5), sharex=False)
    
    for idx, mod in enumerate(["m6A", "m5C", "Y", "I"]):
        ax = axes[idx]
        ax.grid(axis='y', linestyle=':', alpha=0.5, zorder=0)
        
        pos = mapped_positions[mod]
        colors = mod_colors.get(mod, {'line': '#455A64', 'fill': '#CFD8DC'})
        
        if len(pos) > 0:
            try:
                kde = gaussian_kde(pos)
                x_eval = np.linspace(0, 3, 300)
                y_eval = kde(x_eval)
            except:
                counts, edges = np.histogram(pos, bins=100, range=(0, 3), density=True)
                x_eval = (edges[:-1] + edges[1:]) / 2
                y_eval = counts
                
            ax.plot(x_eval, y_eval, color=colors['line'], lw=2.0, zorder=3)
            ax.fill_between(x_eval, 0, y_eval, color=colors['fill'], alpha=0.2, zorder=2)
        
        ax.axvline(1.0, color='#37474F', linestyle='--', lw=1.0, zorder=4)
        ax.axvline(2.0, color='#37474F', linestyle='--', lw=1.0, zorder=4)
        
        ax.set_xlim(0, 3)
        ax.set_ylim(bottom=0)
        
        # Big bold mod label inside subplot
        ax.text(0.05, 0.88, mod, transform=ax.transAxes, fontsize=13, fontweight='bold',
                verticalalignment='top', color=colors['line'],
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor=colors['line'],
                          alpha=0.85, linewidth=1.2), zorder=5)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#37474F')
        ax.spines['bottom'].set_color('#37474F')
        ax.tick_params(colors='#37474F', labelsize=9)
        
    # Set xticks on the bottom subplot only
    axes[-1].set_xticks([0.5, 1.5, 2.5])
    axes[-1].set_xticklabels(["5' UTR", "CDS", "3' UTR"], fontsize=12, fontweight='bold')
    axes[-1].set_xlabel("Transcript Region", fontsize=12, fontweight='bold', labelpad=6)
    
    # Remove x-tick labels from upper panels
    for ax in axes[:-1]:
        ax.set_xticks([0.5, 1.5, 2.5])
        ax.set_xticklabels([])
    
    # Shared y-axis label
    fig.text(0.02, 0.5, "Density", ha="center", va="center", rotation="vertical", fontsize=12, fontweight='bold')
    
    plt.tight_layout(rect=[0.07, 0.08, 0.98, 0.98], h_pad=0.6)
    grid_png = os.path.join(outdir, "panel_f_composite_metagene_grid.png")
    grid_pdf = os.path.join(outdir, "panel_f_composite_metagene_grid.pdf")
    plt.savefig(grid_png, dpi=300)
    plt.savefig(grid_pdf)
    plt.close()
    print(f"Saved 4x1 stacked grid plot to {grid_png}")

    # ----------------- PLOT 3: 1x3 Horizontal Columns (Panel e) -----------------
    # Dimensions match approx 11.5" x 3.2".
    # Fonts are enlarged significantly (ticks/labels=20, labels/titles=22).
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2), sharex=False, sharey=False)
    
    # Calculate y-limit based on m6A density peak
    m6a_pos = mapped_positions.get("m6A", [])
    m6a_max_y = 1.0
    if len(m6a_pos) > 0:
        try:
            kde_m6a = gaussian_kde(m6a_pos)
            x_eval_m6a = np.linspace(0, 3, 300)
            m6a_max_y = np.max(kde_m6a(x_eval_m6a))
        except:
            counts_m6a, _ = np.histogram(m6a_pos, bins=100, range=(0, 3), density=True)
            m6a_max_y = np.max(counts_m6a)
    y_limit_m6a = m6a_max_y * 1.05

    for idx, mod in enumerate(["m6A", "m5C", "Y"]):
        ax = axes[idx]
        ax.grid(axis='y', linestyle=':', alpha=0.5, zorder=0)
        
        pos = mapped_positions[mod]
        colors = mod_colors.get(mod, {'line': '#455A64', 'fill': '#CFD8DC'})
        
        if len(pos) > 0:
            try:
                kde = gaussian_kde(pos)
                x_eval = np.linspace(0, 3, 300)
                y_eval = kde(x_eval)
            except:
                counts, edges = np.histogram(pos, bins=100, range=(0, 3), density=True)
                x_eval = (edges[:-1] + edges[1:]) / 2
                y_eval = counts
                
            ax.plot(x_eval, y_eval, color=colors['line'], lw=2.0, zorder=3)
            ax.fill_between(x_eval, 0, y_eval, color=colors['fill'], alpha=0.2, zorder=2)
        
        ax.axvline(1.0, color='#37474F', linestyle='--', lw=1.0, zorder=4)
        ax.axvline(2.0, color='#37474F', linestyle='--', lw=1.0, zorder=4)
        
        ax.set_xlim(0, 3)
        ax.set_ylim(0, y_limit_m6a)
        
        # Make the mod type name the title of the subplot (fontsize=22, bold, black)
        ax.set_title(mod, fontsize=22, fontweight='bold', pad=10, color='black')
        
        ax.set_xticks([0.5, 1.5, 2.5])
        ax.set_xticklabels(["5' UTR", "CDS", "3' UTR"], fontsize=20, fontweight='bold')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#37474F')
        ax.spines['bottom'].set_color('#37474F')
        ax.tick_params(colors='#37474F', labelsize=20)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
        
    # Shared labels (fontsize=22 and bold)
    fig.text(0.5, 0.02, "Transcript Region", ha="center", va="center", fontsize=22, fontweight='bold')
    fig.text(0.02, 0.5, "Density", ha="center", va="center", rotation="vertical", fontsize=22, fontweight='bold')
    
    plt.tight_layout(rect=[0.05, 0.08, 0.98, 0.98])
    h3col_png = os.path.join(outdir, "panel_f_composite_metagene_horizontal_3col.png")
    h3col_pdf = os.path.join(outdir, "panel_f_composite_metagene_horizontal_3col.pdf")
    plt.savefig(h3col_png, dpi=300)
    plt.savefig(h3col_pdf)
    plt.close()
    print(f"Saved 1x3 horizontal 3col plot to {h3col_png}")

if __name__ == "__main__":
    main()
