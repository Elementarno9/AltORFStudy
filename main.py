import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import plotly.graph_objects as go
import plotly.express as px
from pqdm.processes import pqdm
from tqdm import tqdm
import argparse
import os
import requests

CONFIG = {
    'EMPTY_TAG': 'EMPTY',
    'IMPACTS': {
        'MODIFIER': 0,
        'LOW': 1,
        'MODERATE': 2,
        'HIGH': 3
    }
}


def read_tsv_data(path):
    """Read .tsv file as DataFrame"""
    df = pd.read_csv(path, dtype=str, sep='\t').fillna(CONFIG['EMPTY_TAG'])
    df = df[~(df['ANN[*].IMPACT'] == 'ANN[*].IMPACT')]  # Drop headers
    return df


def read_vcf_data(path):
    """Read .vcf file as DataFrame"""
    df = pd.read_csv(path, sep='\t', header=None).fillna(CONFIG['EMPTY_TAG'])
    return df


def preparing_dataframe(data, drop_duplicates=True):
    """Convert list of dicts into DataFrame"""
    df = pd.DataFrame(data).replace(CONFIG['EMPTY_TAG'], '')
    df['CHROM'] = df['CHROM'].astype('int')
    df['POS'] = df['POS'].astype('int')
    if drop_duplicates:
        df.drop_duplicates(['ID', 'CHROM', 'POS', 'ALT', 'REF', 'ALT_HGVS_C'], inplace=True)
    df = df.sort_values(['CHROM', 'POS']).reset_index(drop=True)

    return df


def save_df(df, path):
    """Rewrite DataFrame to a file"""
    if os.path.exists(path):
        os.remove(path)
    df.to_csv(path, index=False)


def str2list(s):
    """Convert a list-like string into an actual list"""
    for ch in ["'", "[", "]", " "]:
        s = s.replace(ch, '')
    return s.split(',')


def get_rsid(chrom, pos, input_df):
    """Recieve rsID from raw data"""
    res = input_df[(input_df[0] == int(chrom)) & (input_df[1] == int(pos))]
    if len(res) == 0:
        return None
    return res.iloc[0, 2]


def group_snps_by_transcript(df, input_df):
    """Group all rows with SNPs by chromosome & position and divide them by reference and alternative protein"""
    snp_groups = {}

    for i, row in tqdm(df.iterrows(), total=len(df)):
        is_alt = False
        feature_id = row['ANN[*].FEATUREID']
        if feature_id.find('@') != -1:  # refORF
            sep_pos = feature_id.find('@')
        elif feature_id.find('^') != -1:  # alfORF
            sep_pos = feature_id.find('^')
            is_alt = True
        else:  # refORF
            sep_pos = len(feature_id)

        transcript_part = feature_id[:sep_pos]
        protein_part = feature_id[sep_pos + 1:]

        dotpos = transcript_part.find('.')
        if dotpos == -1:
            dotpos = len(transcript_part)
        transcript_id = transcript_part[:dotpos]

        if len(protein_part) == 0:
            protein_part = CONFIG['EMPTY_TAG']

        t = 'alt' if is_alt else 'ref'
        chrom_pos = f'{row["CHROM"]}_{row["POS"]}'

        if chrom_pos not in snp_groups:
            snp_groups[chrom_pos] = {}

        if transcript_id not in snp_groups[chrom_pos]:
            snp_groups[chrom_pos][transcript_id] = {
                'ref': [],
                'alt': []
            }

        snp_groups[chrom_pos][transcript_id][t].append({
            'ID': get_rsid(row['CHROM'], row['POS'], input_df),
            'CHROM': row['CHROM'],
            'POS': row['POS'],
            'REF': row['REF'],
            'ALT': row['ALT'],
            'EFFECT': ' & '.join(str2list(row['ANN[*].EFFECT'])),
            'IMPACT': CONFIG['IMPACTS'][row['ANN[*].IMPACT']],
            'HGVS_P': row['ANN[*].HGVS_P'],
            'HGVS_C': row['ANN[*].HGVS_C'],
            'GENE': row['ANN[*].GENE'],
            'PROTEIN_ACC': protein_part
        })

    return snp_groups


def get_max_impact_snps_iter(args):
    """An iteration of get_max_impact_snps"""
    try:
        chrom_pos, snp_data = args
        max_impacts = []
        for transcript, data in snp_data.items():
            if len(data['alt']) == 0 or len(data['ref']) == 0:
                continue

            ref = data['ref'][0]  # Pick the first
            if len(data['ref']) > 1:
                print(f'WARNING: SNP at chr{ref["CHROM"]}:{ref["POS"]} has {len(data["ref"])} ref proteins.')

            alt_impacts = np.array(list(map(lambda v: v['IMPACT'], data['alt'])))
            max_alt_idxs = np.argwhere(alt_impacts >= 2).reshape(-1)

            if len(max_alt_idxs) == 0:
                max_alt_idxs = np.argmax(alt_impacts).reshape(-1)

            for idx in max_alt_idxs:
                alt = data['alt'][idx]  # Alt with max impact
                max_impacts.append({
                    'ID': ref['ID'],
                    'CHROM': ref['CHROM'],
                    'POS': ref['POS'],
                    'REF': ref['REF'],
                    'ALT': ref['ALT'],
                    'TRANSCRIPT_ACC': transcript,
                    'REF_PROTEIN_ACC': ref['PROTEIN_ACC'],
                    'ALT_PROTEIN_ACC': alt['PROTEIN_ACC'],
                    'REF_IMPACT': ref['IMPACT'],
                    'ALT_IMPACT': alt['IMPACT'],
                    'REF_HGVS_P': ref['HGVS_P'],
                    'ALT_HGVS_P': alt['HGVS_P'],
                    'REF_HGVS_C': ref['HGVS_C'],
                    'ALT_HGVS_C': alt['HGVS_C'],
                    'GENE': ref['GENE'],
                    'REF_EFFECT': ref['EFFECT'],
                    'ALT_EFFECT': alt['EFFECT'],
                })
        return max_impacts
    except Exception as e:
        print(f'ERROR IN LOOP: {e}')
        return []


def get_max_impact_snps(snp_groups, parallel=True):
    """Choosing from each SNP group only with max impact (>= 2 or at least anyone)"""
    if parallel:
        n_jobs = max(1, os.cpu_count() - 1)
        print(f'Paralleling into {n_jobs} processes...')
        results = pqdm(snp_groups.items(), get_max_impact_snps_iter, n_jobs=n_jobs)
    else:
        results = [get_max_impact_snps_iter(args) for args in tqdm(snp_groups.items(), total=len(snp_groups))]

    return [res for sublist in results for res in sublist]


def add_gtex_info(df):
    """Make request to gtexportal.org to find associations for SNPs"""
    tissue_col, nes_col, pvalue_col, idx = [], [], [], []

    for i, row in tqdm(df.iterrows(), total=len(df)):
        rsid = row['ID']
        req = requests.get(
            (f'https://gtexportal.org/rest/v1/association/singleTissueEqtl?format=json&snpId={rsid}&'
             'tissueSiteDetailId=Brain_Amygdala%2CBrain_Anterior_cingulate_cortex_BA24%2CBrain_Caudate_basal_ganglia%2C'
             'Brain_Cerebellar_Hemisphere%2CBrain_Cerebellum%2CBrain_Cortex%2CBrain_Frontal_Cortex_BA9%2C'
             'Brain_Hippocampus%2CBrain_Hypothalamus%2CBrain_Nucleus_accumbens_basal_ganglia%2C'
             'Brain_Putamen_basal_ganglia%2CBrain_Spinal_cord_cervical_c-1%2C'
             'Brain_Substantia_nigra&datasetId=gtex_v8')).json()

        if len(req['singleTissueEqtl']) == 0:
            continue
        best_tissue = min(req['singleTissueEqtl'], key=lambda x: x['pValue'])
        tissue_col.append(best_tissue['tissueSiteDetailId'])
        nes_col.append(best_tissue['nes'])
        pvalue_col.append(best_tissue['pValue'])
        idx.append(i)

    tissue_df = pd.DataFrame({'Tissue': tissue_col, 'NES': nes_col, 'Tissue p-value': pvalue_col}, index=idx)
    return df.join(tissue_df)


def plot_sankey_of_effects(df, base_path, tag):
    """Plot sankey of changing of SNPs' effects"""
    effect_grouped = df.groupby(['REF_EFFECT', 'ALT_EFFECT']).size().reset_index()
    ref_count, alt_count = [effect_grouped.groupby(e).sum().reset_index() for e in ['REF_EFFECT', 'ALT_EFFECT']]
    effect_unique = pd.unique(pd.concat([effect_grouped['REF_EFFECT'], effect_grouped['ALT_EFFECT']]))

    le = LabelEncoder()
    le.fit(effect_unique)
    source = le.transform(effect_grouped['REF_EFFECT'])
    target = le.transform(effect_grouped['ALT_EFFECT']) + len(le.classes_)
    weights = effect_grouped[0]
    labels = le.classes_.tolist() * 2

    for i, label in enumerate(labels):
        if i < len(le.classes_):
            row = ref_count[ref_count['REF_EFFECT'] == label]
        else:
            row = alt_count[alt_count['ALT_EFFECT'] == label]

        count = int(row[0]) if len(row) == 1 else 0

        labels[i] = f'{label} ({count})'

    node_colors = px.colors.qualitative.Set1 + px.colors.qualitative.Set2 + px.colors.qualitative.Set3
    link_colors = np.array((px.colors.qualitative.Pastel1 + px.colors.qualitative.Pastel2) * 2)
    link_colors = link_colors[source]

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20,
            thickness=15,
            line=dict(color="black", width=0.5),
            label=labels,
            color=node_colors
        ),
        link=dict(
            source=source,
            target=target,
            value=weights,
            color=link_colors
        ))])

    fig.update_traces(name='label+value')

    fig.update_layout(annotations=[{
        'xref': 'paper',
        'yref': 'paper',
        'x': 0,
        'y': 1.1,
        'text': 'Effect in refORFs',
        'showarrow': False,
        'font': {'size': 20, 'color': 'black'}
    }, {
        'xref': 'paper',
        'yref': 'paper',
        'x': 1,
        'y': 1.1,
        'text': 'Effect in altORFs',
        'showarrow': False,
        'font': {'size': 20, 'color': 'black'}
    }
    ])
    fig.update_layout(title_text="Changing effect between reference & alternative ORFs",
                      title_font_size=32,
                      font_size=14)

    fig.write_image(os.path.join(base_path, f"sankey_effect_{tag}.png"), width=1000, height=600)
    fig.write_html(os.path.join(base_path, f"sankey_effect_{tag}.html"))
    fig.show()


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--input', required=True, help="Path to the folder with OpenVar data.")
    parser.add_argument('--output', help="Path to the output folder. Default is same as input path.")
    parser.add_argument('--name', help="The prefix for .tsv file, if it is not the same as input path.")
    parser.add_argument('--parallel', action='store_true', help="Using a few cores for faster run.")
    parser.add_argument('--only-sankey', action='store_true', help=("Only replot sankey for existing output."
                                                                    "Can be used when plotting parameters need to"
                                                                    "be changed."))

    args = parser.parse_args()
    if args.name is None:
        args.name = args.input
    if args.output is None:
        args.output = args.input

    # Load data
    INPUT_PATH = os.path.join(args.input, 'output', 'input_vcf.vcf')
    DATA_PATH = os.path.join(args.input, 'output', f'{args.name}_annOnePerLine.tsv')
    OUTPUT_PATH = args.output

    if args.only_sankey:
        filename = os.path.join(OUTPUT_PATH, 'top_impact.csv')
        print(f'Plotting sankey for {filename}')
        top_impact_df = pd.read_csv(filename).fillna('')
        plot_sankey_of_effects(top_impact_df, OUTPUT_PATH, 'top')
        return

    if not os.path.exists(DATA_PATH):
        print(f"ERROR: incorrent name `{args.name}` in path to .tsv file.")
        parser.print_help()
        exit()

    df = read_tsv_data(DATA_PATH)
    input_df = read_vcf_data(INPUT_PATH)

    print(f'Data loaded from {DATA_PATH}. Total {len(df)} rows.')

    # Main part
    snp_groups = group_snps_by_transcript(df, input_df)
    print(f'Grouped into {len(snp_groups)} different variants.')
    max_impact_snps = get_max_impact_snps(snp_groups, parallel=args.parallel)

    if len(max_impact_snps) == 0:
        print("ERROR: no SNPs in final result. Check is data available.")
        exit()

    # Prepare output
    print(f'Preparing output of {len(max_impact_snps)} SNPs.')
    impact_df = preparing_dataframe(max_impact_snps)
    top_impact_df = impact_df[impact_df['ALT_IMPACT'] >= 2]

    print("Adding information from GTEx...")
    top_impact_df = add_gtex_info(top_impact_df)

    print(f"Saving data.")
    save_df(impact_df, os.path.join(OUTPUT_PATH, 'impact.csv'))
    save_df(top_impact_df, os.path.join(OUTPUT_PATH, 'top_impact.csv'))

    # Plot graphics
    if len(top_impact_df) > 0:
        plot_sankey_of_effects(top_impact_df, OUTPUT_PATH, 'top')


if __name__ == '__main__':
    main()
