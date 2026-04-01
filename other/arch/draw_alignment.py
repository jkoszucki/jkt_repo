from pathlib import Path
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, ConnectionPatch
from matplotlib.colors import to_rgb, to_hex

# --- UTILITY FUNCTIONS ---

def darken_color(hex_color, factor=0.7):
    """
    Darken a hex color by a given factor (0 < factor < 1).
    Returns a new hex color string.
    """
    rgb = np.array(to_rgb(hex_color))
    dark_rgb = (rgb * factor).clip(0, 1)
    return to_hex(dark_rgb)


def load_category_map(path):
    """
    Load a two-column TSV mapping proteinID to category/source (e.g. 'GENSCRIPT', 'PREDICTION', etc.).
    Returns { proteinID: CATEGORY }.
    """
    category_map = {}
    with open(path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            protein, source = parts[0].strip(), parts[1].strip().upper()
            category_map[protein] = source
    return category_map


# --- COLOR DICT GENERATOR ---

def generate_subject_color_dict(
    predictions_and_enzymes_path,
    protein_ids,
    category_colors,
    default_color='salmon'
):
    """
    Generate mapping from proteinID to hex color based on prediction_strength or category.
    Priority:
      1. prediction_strength 'strong' -> #86c187
      2. prediction_strength 'likely' -> #dbecdc
      3. category_colors based on category map (e.g. 'GENSCRIPT' => blue)
      4. default_color
    """
    df_full = pd.read_csv(predictions_and_enzymes_path, sep='\t', dtype=str)
    category_map = load_category_map(predictions_and_enzymes_path)
    color_map = {}
    for pid in protein_ids:
        rows = df_full[df_full['proteinID'] == pid]
        for _, row in rows.iterrows():
            strength = str(row.get('prediction_strength', '')).strip().lower()
            if strength == 'strong':
                color = '#86c187'
            elif strength == 'likely':
                color = '#dbecdc'
            else:
                cat = category_map.get(pid, '')
                color = category_colors.get(cat, default_color)
            color_map[pid] = color
            break
    return color_map


# --- PLOTTING FUNCTION WITH HSP MERGING & REORDERING ---

def plot_query_vs_reference(
    predictions_and_enzymes,
    protein_ids,
    query_id,
    workdir,
    cc_key,
    tmp_folder=None,
    fname=None,
    row_spacing=1.0,
    block_height=0.4,
    edge_linewidth=1,
    connector_linewidth=1.5,
    evalue_threshold=1e-3,
    subject_colors=None,
    hsp_merge_gap=50
):
    """
    Draw one query bar at the top, then for each other protein (subject):

    1) Merge HSPs whose query‐side intervals are within `hsp_merge_gap` (default 50 aa),
       computing a weighted‐average percent identity.
    2) Reorder subjects so that:
         (a) All "blue" proteins (category == 'GENSCRIPT') appear first,
         (b) Then those sharing the same KL as the query (sorted by descending average pident),
         (c) Then the rest, grouped by KL (ascending), each group sorted by descending pident.
    3) Plot the merged aligned region (plus flanks) under the query, 
       and draw a curved connector from the query bar to the subject bar, 
       labeling it with the merged percent identity.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    tmp_folder = Path(tmp_folder).expanduser() if tmp_folder else Path.home() / 'Downloads' / 'blastp_tmp'
    tmp_folder.mkdir(parents=True, exist_ok=True)

    # Load sequences + metadata
    if isinstance(predictions_and_enzymes, pd.DataFrame):
        df = predictions_and_enzymes.copy()
    else:
        df = pd.read_csv(predictions_and_enzymes, sep='\t', dtype=str)

    df = df[df['proteinID'].isin(protein_ids)].reset_index(drop=True)
    if query_id not in df['proteinID'].tolist():
        raise ValueError("query_id must be in protein_ids")

    # Build a component-specific assign_map so that duplicate proteinIDs
    # (e.g. PC0692) get the correct KL per connected component:
    assign_map = {}
    for pid in protein_ids:
        # Subset rows for this pid
        rows = df[df['proteinID'] == pid]
        if len(rows) > 1 and pid == 'PC0692':
            # Two entries for PC0692: KL111 for cc4, KL127 for cc14
            if cc_key == 'cc4':
                assign_map[pid] = 'KL111'
            elif cc_key == 'cc14':
                assign_map[pid] = 'KL127'
            else:
                # If somehow appears elsewhere, pick the first available KL
                assign_map[pid] = rows.iloc[0]['assing_K_locus_host_when_specificity_missing']
        else:
            # Only one row or not PC0692: take the KL directly
            assign_map[pid] = rows.iloc[0]['assing_K_locus_host_when_specificity_missing']

    # category_map: maps proteinID -> source (e.g. 'GENSCRIPT', 'PREDICTION', etc.)
    category_map = load_category_map(predictions_and_enzymes)
    seqs = dict(zip(df['proteinID'], df['seq']))
    lengths = {pid: len(seq) for pid, seq in seqs.items()}

    # Write FASTA for BLASTp
    fasta = tmp_folder / 'sequences.fasta'
    with open(fasta, 'w') as fh:
        for _, row in df.iterrows():
            fh.write(f">{row['proteinID']}\n{row['seq']}\n")

    blast_out = tmp_folder / 'blastp_raw.tsv'
    subprocess.run(['makeblastdb', '-in', str(fasta), '-dbtype', 'prot'], check=True)
    subprocess.run([
        'blastp', '-query', str(fasta), '-db', str(fasta),
        '-out', str(blast_out), '-outfmt',
        '6 qseqid sseqid evalue pident bitscore qlen qstart qend slen sstart send'
    ], check=True)

    # Parse HSPs
    cols = ['query','target','evalue','pident','bitscore',
            'qlen','qstart','qend','slen','sstart','send']
    edges = pd.read_csv(blast_out, sep='\t', names=cols, header=None)
    edges = edges[edges['evalue'] <= evalue_threshold]

    # Build a DataFrame of hits where one side is query_id
    hit_list = []
    for _, row in edges.iterrows():
        if row['query'] == query_id and row['target'] in protein_ids and row['target'] != query_id:
            hit_list.append({
                'subj': row['target'],
                'qstart': int(row['qstart']),
                'qend':   int(row['qend']),
                'sstart': int(row['sstart']),
                'send':   int(row['send']),
                'pident': float(row['pident']),
                'bitscore': float(row['bitscore'])
            })
        elif row['target'] == query_id and row['query'] in protein_ids and row['query'] != query_id:
            hit_list.append({
                'subj': row['query'],
                'qstart': int(row['sstart']),
                'qend':   int(row['send']),
                'sstart': int(row['qstart']),
                'send':   int(row['qend']),
                'pident': float(row['pident']),
                'bitscore': float(row['bitscore'])
            })
    hits_df = pd.DataFrame(hit_list)
    if hits_df.empty:
        raise ValueError(f"No BLASTp hits involving {query_id} at evalue<={evalue_threshold}")

    # Merge HSPs for each subject if their query‐side gap < hsp_merge_gap
    merged_rows = []
    for subj, group in hits_df.groupby('subj'):
        sub_df = group.sort_values('bitscore', ascending=False).reset_index(drop=True)
        top = sub_df.loc[0].copy()
        merged = top.copy()

        for idx in range(1, len(sub_df)):
            other = sub_df.loc[idx]
            # Check non-overlap on query side
            if other['qend'] < merged['qstart'] or other['qstart'] > merged['qend']:
                # Compute gap on query side
                if other['qstart'] > merged['qend']:
                    gap = other['qstart'] - merged['qend'] - 1
                else:
                    gap = merged['qstart'] - other['qend'] - 1

                if gap < hsp_merge_gap:
                    # Merge these two HSPs
                    new_qstart = min(merged['qstart'], other['qstart'])
                    new_qend   = max(merged['qend'],   other['qend'])
                    len1 = merged['qend'] - merged['qstart'] + 1
                    len2 = other['qend']   - other['qstart']   + 1
                    avg_pident = (merged['pident'] * len1 + other['pident'] * len2) / (len1 + len2)
                    new_sstart = min(merged['sstart'], other['sstart'])
                    new_send   = max(merged['send'],   other['send'])
                    merged['qstart'] = new_qstart
                    merged['qend']   = new_qend
                    merged['pident'] = avg_pident
                    merged['sstart'] = new_sstart
                    merged['send']   = new_send
                    break

        merged_rows.append({
            'subj':    subj,
            'qstart':  merged['qstart'],
            'qend':    merged['qend'],
            'sstart':  merged['sstart'],
            'send':    merged['send'],
            'pident':  merged['pident']
        })

    merged_df = pd.DataFrame(merged_rows)

    # Add columns: KL, same_kl, category, is_blue
    merged_df['KL'] = merged_df['subj'].map(assign_map)
    query_KL = assign_map.get(query_id, '')
    merged_df['same_kl'] = merged_df['KL'] == query_KL

    merged_df['category'] = merged_df['subj'].map(category_map).fillna('')
    # “Blue proteins” are those with category == 'GENSCRIPT'
    merged_df['is_blue'] = merged_df['category'] == 'GENSCRIPT'

    # Sort order:
    # 1. is_blue=True first
    # 2. same_kl=True next
    # 3. then by KL (ascending) for the remainder
    # 4. within each group, by descending pident
    merged_df_sorted = merged_df.sort_values(
        by=['is_blue', 'same_kl', 'KL', 'pident'],
        ascending=[False, False, True, False]
    ).reset_index(drop=True)

    partners = merged_df_sorted['subj'].tolist()
    n = len(partners)

    # Prepare figure
    total_height = (n + 1) * row_spacing + row_spacing
    fig, ax = plt.subplots(figsize=(8, max(3, total_height)), dpi=300)

    query_len = lengths[query_id]
    query_color = subject_colors.get(query_id, '#86c187') if subject_colors else '#86c187'
    query_edge = darken_color(query_color)

    # Y position for the single query bar
    y_query = (n + 1) * row_spacing

    # 1) Draw the single query bar
    ax.add_patch(Rectangle((0, y_query - block_height / 2), query_len, block_height,
                           facecolor=query_color, edgecolor=query_edge,
                           linewidth=edge_linewidth, zorder=2))

    # Prepare x‐positions for connectors, evenly spaced
    x_margin = max(1, int(0.02 * query_len))
    if n > 1:
        x_connectors = np.linspace(x_margin, query_len - x_margin, n)
    else:
        x_connectors = np.array([query_len / 2.0])

    # Collect Y‐ticks & labels
    yticks = [y_query]
    ylabels = [f"{query_id} ({assign_map.get(query_id,'')})"]

    # Draw each subject row & connector in sorted order
    for i, subj in enumerate(partners):
        row = merged_df_sorted.iloc[i]
        y_sub = (n - i) * row_spacing
        yticks.append(y_sub)
        ylabels.append(f"{subj} ({assign_map.get(subj,'')})")

        subj_len = lengths[subj]
        subj_color = subject_colors.get(subj, 'salmon') if subject_colors else 'salmon'
        subj_edge = darken_color(subj_color)

        q_start = int(row['qstart'])
        q_end   = int(row['qend'])
        s_start = int(row['sstart'])
        s_end   = int(row['send'])
        pident  = float(row['pident'])

        # Clip query‐side coordinates
        q_start = max(1, q_start)
        q_end = min(query_len, q_end)

        # Compute subject flanking lengths
        left_len = s_start - 1
        right_len = subj_len - s_end

        # LEFT FLANK under query
        max_left_fit = q_start - 1
        left_width = min(left_len, max_left_fit)
        left_x = (q_start - 1) - left_width

        # MIDDLE ALIGNED BLOCK
        mid_x = q_start - 1
        mid_width = q_end - q_start + 1

        # RIGHT FLANK under query
        max_right_fit = query_len - q_end
        right_width = min(right_len, max_right_fit)
        right_x = q_end

        # Draw left flank if any
        if left_width > 0:
            ax.add_patch(Rectangle((left_x, y_sub - block_height / 2),
                                   left_width, block_height,
                                   facecolor='lightgray', alpha=0.4,
                                   edgecolor=None, zorder=1))

        # Draw aligned subject block
        ax.add_patch(Rectangle((mid_x, y_sub - block_height / 2),
                               mid_width, block_height,
                               facecolor=subj_color, edgecolor=subj_edge,
                               linewidth=edge_linewidth, zorder=2))

        # Draw right flank if any
        if right_width > 0:
            ax.add_patch(Rectangle((right_x, y_sub - block_height / 2),
                                   right_width, block_height,
                                   facecolor='lightgray', alpha=0.4,
                                   edgecolor=None, zorder=1))

        # Draw curved connector on top
        x_conn = x_connectors[i]
        y_conn_query = y_query - block_height / 2
        y_conn_sub = y_sub + block_height / 2
        midpoint = query_len / 2.0
        rad = 0.3 if x_conn < midpoint else -0.3  # swapped curvature

        ax.add_patch(ConnectionPatch(
            (x_conn, y_conn_query),
            (x_conn, y_conn_sub),
            coordsA='data', coordsB='data',
            arrowstyle='-',
            linewidth=connector_linewidth,
            connectionstyle=f'arc3,rad={rad}',
            color='gray', alpha=0.3,
            zorder=5
        ))

        # Label merged percent identity next to connector end
        ax.text(
            x_conn + 50, y_conn_sub + 0.02,
            f"{pident:.1f}%",
            fontweight='bold', ha='center', va='bottom',
            fontsize=12, color='black', zorder=6
        )

    # Set y‐axis ticks & labels
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontweight='bold', fontsize=8)
    ax.invert_yaxis()

    # X‐axis margin
    ax.set_xlim(-x_margin, query_len + x_margin)

    # Y‐limits with margin
    ax.set_ylim(0.5, (n + 1) * row_spacing + row_spacing - 0.5)

    # Style axes
    ax.set_xlabel('Residue Position on Query')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Adjust left margin (so y‐labels are not truncated)
    fig.subplots_adjust(left=0.25, right=0.95, top=0.95, bottom=0.05)

    # Adjust font size and weight for axis labels and ticks
    ax.set_xlabel(ax.get_xlabel(), fontsize=12, fontweight='bold')
    ax.yaxis.set_tick_params(labelsize=12)
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')

    # Save
    out_name = fname if fname else f"{query_id}_alignment_pairs"
    out_path = Path(workdir) / f"{out_name}.pdf"
    fig.tight_layout()
    fig.savefig(out_path, format='pdf', dpi=300)
    plt.close(fig)
    return str(out_path)


# --- EXECUTION BLOCK ---

if __name__ == "__main__":
    predictions_and_enzymes = (
        '/Users/januszkoszucki/MGG Dropbox/Projects/'
        'kleb-prophage-div/2025-02-12_KLEBDATA_LIGHT/'
        '3_GWAS/2025-05-13_KLEBDATA_LIGHT_GWAS/'
        '4_ANALYZE/2_AGGREGATED_DATA/predictions_and_enzymes.tsv'
    )

    workdir = (
        '/Users/januszkoszucki/MGG Dropbox/Projects/'
        'kleb-prophage-div/2025-02-12_KLEBDATA_LIGHT/'
        '3_GWAS/SUPPLEMENT/'
        'supplementary_text_S3/1_ALIGNMENTS/'
    )

    tmp_folder = '/Users/januszkoszucki/Downloads/blastp_tmp_dir'
    category_colors = {
        'GENSCRIPT': '#6baed6',        # blue
        'PREDICTION': '#86c187',
        'LITERATURE_SEARCH': '#b695c5',
        'PROPHAGE_ZDKLAB': '#de69a5'
    }

    connected_components_dict = {
        'cc1': {'protein_ids': ['PC0279','PC0449','1251_37','PC1397','391_03','PC0406'], 'query_id': 'PC0449'},
        'cc2': {'protein_ids': ['PC0772','1091_44','PC0871','1723_59','1724_71','BBF66868_1','YP_009153197_1'], 'query_id': 'PC0772'},
        'cc3': {'protein_ids': ['YP_009153202_1','PC1103','PC1375','0574_17','QOV05454_1','YP_009797016_1'], 'query_id': 'PC1103'},
        'cc4': {'protein_ids': ['PC1605','YP_009153199_1','PC1712','PC0692','PC0538','184_04'], 'query_id': 'PC1605'},
        'cc5': {'protein_ids': ['0496_72','YP_009153200_1','PC1357', '248_38'], 'query_id': 'PC1357'},
        'cc6': {'protein_ids': ['PC0844','PC0551','1863_75','1630_70','145_08'], 'query_id': 'PC0844'},
        'cc7': {'protein_ids': ['021_14','PC0367'], 'query_id': 'PC0367'},
        'cc8': {'protein_ids': ['YP_003347555_1','UCR74082_1','PC0537'], 'query_id': 'PC0537'},
        'cc12': {'protein_ids': ['PC0827','APZ82804_1','YP_009153203_1'], 'query_id': 'PC0827'},
        'cc13': {'protein_ids': ['0367_12','PC0619','WOK01638_1','0391_11','PC0671','914_74'], 'query_id': 'PC0619'},
        'cc14': {'protein_ids': ['319_37','PC0959','PC0692','PC1803'], 'query_id': 'PC0692'},
        'cc15': {'protein_ids': ['UCR74083_1','617_77','PC0854','BAN78446_1'], 'query_id': 'PC0854'},
        'cc16': {'protein_ids': ['1409_59', '1248_57', '1441_47', 'QOI68577_1', 'QIW88225_1'], 'query_id': '1409_59'}
    }

    for fname, data in connected_components_dict.items():
        protein_ids = data['protein_ids']
        query_id = data['query_id']
        subject_colors = generate_subject_color_dict(
            predictions_and_enzymes, protein_ids, category_colors, default_color='salmon'
        )
        out_path = plot_query_vs_reference(
            predictions_and_enzymes,
            protein_ids,
            query_id,
            workdir,
            cc_key=fname,
            tmp_folder=tmp_folder,
            fname=fname,
            row_spacing=0.7,
            block_height=0.3,
            connector_linewidth=1.5,
            evalue_threshold=1e-3,
            subject_colors=subject_colors,
            hsp_merge_gap=100
        )
        print("Plot saved to:", out_path)
