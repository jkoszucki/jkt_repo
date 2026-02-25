import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the metadata
df = pd.read_csv('/Users/januszkoszucki/MGG Dropbox/Projects/kleb-prophage-div/2025-02-12_KLEBDATA_LIGHT/1_BACTERIA/bacteria_metadata.tsv', sep='\t')

def assign_category(row):
    """Assign isolates to 'Disease', 'Clinical', or 'Environmental'."""
    if str(row.get('kaspah_infection','')).lower() == 'yes':
        return 'Disease'
    ks = str(row.get('klebpavia_group_summary','')).lower()
    if 'human - hospital disease' in ks:
        return 'Disease'
    for col in ['kaspah_infection','kaspah_specimen']:
        if pd.notnull(row.get(col)) and row[col] not in ['', 'UNK']:
            return 'Clinical'
    if 'human - hospital' in ks:
        return 'Clinical'
    return 'Environmental'

def plot_KSC_diagonal_custom(df, species_abbr_list, min_sc_count=20, min_k_count=20,
                             jitter=0.3, output_file='KSC.pdf'):
    """
    Custom diagonal sort:
      1. Exclude UNK for sorting
      2. SCs with one K: sort by isolate count (highest last)
      3. Place their single K on the diagonal
      4. SCs with two Ks: place the more frequent K on diagonal (or the other if already placed)
      5. SCs with >2: place their top Ks on next free diagonal slots
    """
    df_sub = df[df['species_abbreviation'].isin(species_abbr_list)].copy()
    df_sub['category'] = df_sub.apply(assign_category, axis=1)

    # Filter SCs and K loci by count
    sc_counts = df_sub['MGG_SC'].value_counts()
    keep_scs = sc_counts[sc_counts>=min_sc_count].index.tolist()
    df_sc = df_sub[df_sub['MGG_SC'].isin(keep_scs)]
    k_counts = df_sc['MGG_K_locus'].value_counts()
    keep_ks_all = k_counts[k_counts>=min_k_count].index.tolist()

    # Exclude UNK for sorting
    keep_ks = [k for k in keep_ks_all if k!='UNK']
    keep_scs = [s for s in keep_scs if s!='UNK']

    # contingency
    pc = df_sc[df_sc['MGG_K_locus'].isin(keep_ks)].groupby(['MGG_SC','MGG_K_locus']).size().unstack(fill_value=0)

    # unique counts
    uniqK = (pc>0).sum(axis=1)
    # isolate counts per cell in pc

    sc_single = [sc for sc in uniqK.index if uniqK[sc]==1]
    sc_double = [sc for sc in uniqK.index if uniqK[sc]==2]
    sc_multi  = [sc for sc in uniqK.index if uniqK[sc]>2]

    # sort singles by total isolates in that SC (descending at bottom -> ascending)
    total_sc = pc.sum(axis=1)
    sc_single_sorted = sorted(sc_single, key=lambda sc: total_sc[sc])

    # build diagonal K order
    k_order = []
    sc_order = []

    # step 2 & 3: singles
    for sc in sc_single_sorted:
        sc_order.append(sc)
        # find its single K
        k = pc.loc[sc][pc.loc[sc]>0].idxmax()
        if k not in k_order:
            k_order.append(k)

    # step 4: doubles
    for sc in sorted(sc_double, key=lambda sc: total_sc[sc]):
        sc_order.append(sc)
        row = pc.loc[sc]
        # sort its two Ks by count desc
        two = row[row>0].sort_values(ascending=False)
        # place primary
        k1, k2 = list(two.index[:2])
        if k1 not in k_order:
            k_order.append(k1)
        elif k2 not in k_order:
            k_order.append(k2)

    # step 5: multis (e.g. top two per SC)
    for sc in sorted(sc_multi, key=lambda sc: total_sc[sc]):
        sc_order.append(sc)
        row = pc.loc[sc].sort_values(ascending=False)
        for k in row.index[:2]:
            if k not in k_order:
                k_order.append(k)

    # append any remaining
    for k in keep_ks_all:
        if k not in k_order:
            k_order.append(k)
    # add UNK last if present
    if 'UNK' in keep_ks_all:
        k_order.append('UNK')
    if 'UNK' in df_sc['MGG_SC'].unique():
        sc_order.append('UNK')

    # prepare plotting df
    df_plot = df_sc[df_sc['MGG_K_locus'].isin(k_order)].copy()
    df_plot['sc_idx'] = df_plot['MGG_SC'].map({s:i for i,s in enumerate(sc_order)})
    df_plot['k_idx'] = df_plot['MGG_K_locus'].map({k:i for i,k in enumerate(k_order)})

    # jitter
    np.random.seed(42)
    df_plot['x_jit'] = df_plot['k_idx'] + np.random.uniform(-jitter, jitter, len(df_plot))
    df_plot['y_jit'] = df_plot['sc_idx'] + np.random.uniform(-jitter, jitter, len(df_plot))

    # plot
    ncols, nrows = len(k_order), len(sc_order)
    fig, ax = plt.subplots(figsize=(max(8,ncols*0.5+3), max(6,nrows*0.4+3)))
    color_map={'Disease':'red','Clinical':'orange','Environmental':'green'}
    for cat,grp in df_plot.groupby('category'):
        ax.scatter(grp['x_jit'], grp['y_jit'], s=5, alpha=0.6, color=color_map[cat], label=cat)

    ax.set_xticks(range(ncols)); ax.set_xticklabels(k_order, rotation=90, va='top')
    ax.set_yticks(range(nrows)); ax.set_yticklabels(sc_order)
    ax.set_xlim(-0.5, ncols-0.5); ax.set_ylim(-0.5, nrows-0.5)
    ax.set_aspect('equal')
    ax.set_xticks(np.arange(-0.5,ncols,1), minor=True)
    ax.set_yticks(np.arange(-0.5,nrows,1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
    ax.set_xlabel('K locus'); ax.set_ylabel('SC')
    ax.set_title('K locus vs SC (Custom Diagonal Sort)')
    ax.legend(title='Category', bbox_to_anchor=(1.05,1), loc='upper left')
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved custom diagonal plot to {output_file}")

# Run for K. pneumoniae complex
plot_KSC_diagonal_custom(df, ['KPN','KVV','KQQ','KQS'], min_sc_count=20, min_k_count=20)
