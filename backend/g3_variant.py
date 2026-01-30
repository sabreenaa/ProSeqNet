# COQ8B Protein Analysis Project
# Tasks 3 & 4: Genetic Variants Mapping and Data Visualization
import os
import pandas as pd
import numpy as np
import requests
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import warnings
warnings.filterwarnings('ignore')

import requests
import textwrap
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns

def fetch_variant_data(uniprot_id):
    uniprot_id = uniprot_id.strip().upper()
    url = f'https://www.ebi.ac.uk/proteins/api/variation/{uniprot_id}?format=json'
    # This is the API link from Uniprot
    #print('Requesting', url)
    # Fetch the JSON
    r = requests.get(url, timeout=30)
    # if 404 or not found, raise a clear error
    if r.status_code == 404:
        raise ValueError(f"No data found for UniProt ID {uniprot_id} (404)")
    r.raise_for_status()
    try:
        data = r.json()
    except ValueError:
        raise ValueError(f"Invalid JSON response for UniProt ID {uniprot_id}")
    variants = data.get('features')
    if not variants:
        raise ValueError(f"No variant data ('features') found for UniProt ID {uniprot_id}")
    return variants

def variant_dataframe(uniprot_id):
    variants = fetch_variant_data(uniprot_id)
    parsed_variants = []
    for var in variants:
        if var.get('type') == 'VARIANT':
            xrefs = var.get('xrefs') or var.get('xref') or []
            first_id = None
            external_link = None
            for x in xrefs:
                if not isinstance(x, dict):
                    continue
                if first_id is None and 'id' in x:
                    first_id = x['id']
                if external_link is None and 'url' in x:
                    external_link = x['url']

            parsed_variants.append({
                'variant_id': first_id,
                'external_url': external_link,  
                'type': var.get('type'),
                'alt_seq': var.get('alternativeSequence'),
                'begin': var.get('begin'),
                'end': var.get('end'),
                'genomicLocation': var.get('genomicLocation'),
                'consequence': var.get('consequenceType'),
                'mutatedType': var.get('mutatedType'),
                'predictions': var.get('predictions'),
                'wild_type': var.get('wildType'),
                'association' : var.get('association'),
                'Clinical Significance': var.get('clinicalSignificance'),
            })
    # Convert to DataFrame
    df_variants = pd.DataFrame(parsed_variants)
    # Predicted Dataframe
    prediction_records = []
    for v in variants:
        if v.get('type') == 'VARIANT':
            position = v.get("begin")
            type = v.get('type')
            alt_seq = v.get('alternativeSequence')
            end = v.get('end')
            genomicLocation = v.get('genomicLocation')
            consequence= v.get('consequenceType')
            mutatedType= v.get('mutatedType')
            wild_type = v.get('wildType')
            xrefs = v.get('xref',[])
            first_id = None
            xrefs = v.get('xrefs', [])
            for x in xrefs:
                if 'id' in x:
                    first_id = x['id']
                break
            for p in v.get('predictions', []):
                prediction_records.append({
                    'variant_id': first_id,
                    'position': position,
                    'type': type,
                    'alt_seq': alt_seq,
                    'end': end,
                    'genomicLocation': genomicLocation,
                    'consequence': consequence,
                    'mutatedType': mutatedType,
                    'wild_type': wild_type,
                    'algorithm': p.get('predAlgorithmNameType'),
                    'prediction': p.get('predictionValType'),
                    'score': p.get('score'),
                    'source': ",".join(p.get('sources', []))
                })

    pred_df = pd.DataFrame(prediction_records)
    # clean up preditions
    impact_map = {
        "deleterious": "Deleterious",
        "deleterious - low confidence": "Deleterious (low conf.)",
        "tolerated": "Tolerated",
        "tolerated - low confidence": "Tolerated (low conf.)",
        "benign": "Benign",
        "possibly damaging": "Possibly Damaging",
        "probably damaging": "Probably Damaging"
    }

    def map_prediction(pred):
        pred = pred.lower()
        if "probably damaging" in pred or "deleterious" in pred:
            return "High impact"
        elif "possibly damaging" in pred:
            return "Moderate impact"
        elif "tolerated" in pred or "benign" in pred:
            return "Low / neutral"
        elif "low confidence" in pred or "unknown" in pred:
            return "Uncertain"
        else:
            return "Other"

    pred_df["impact_class"] = pred_df["prediction"].apply(map_prediction)
    pred_df["Map_pred"] = (
        pred_df["prediction"]
        .str.lower()
        .map(lambda x: next((impact_map[k] for k in impact_map if k in x), "Other"))
    )

    pred_df["position"] = pd.to_numeric(pred_df["position"], errors="coerce")
    pred_df['score'] = pd.to_numeric(pred_df['score'], errors='coerce')

    return df_variants, pred_df

############# PLOTTING FUNCTION #######################

def Variant_analysis(uniprot_id):
    df_variants, pred_df = variant_dataframe(uniprot_id)
    polyphen_df = pred_df[pred_df['algorithm'].astype(str).str.contains('polyphen', case=False, na=False)].copy()
    polyphen_df['score'] = pd.to_numeric(polyphen_df['score'], errors='coerce')

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
    fig = plt.figure(figsize=(16, 10))
    fig
    # A
    plt.subplot(2, 3, 1)
    plt.hist(polyphen_df['position'].dropna(), bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel("Amino acid position")
    plt.ylabel("Number of variants")
    plt.title("A. Variant distribution along protein sequence")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    # B — show source (SIFT vs PolyPhen) by color/marker
    plt.subplot(2, 3, 2)

    BIN_SIZE = 20

    pred_df1 = pred_df.copy()
    pred_df1["pos_bin"] = (pred_df1["position"] // BIN_SIZE) * BIN_SIZE

    agg_df = (
        pred_df1
        .groupby(["pos_bin", "Map_pred", "algorithm"])
        .size()
        .reset_index(name="count")
    )

    impact_order = [
        "Probably Damaging",
        "Possibly Damaging",
        "Benign",
        "Deleterious",
        "Tolerated"
    ]

    y_map = {label: i for i, label in enumerate(impact_order)}
    agg_df["y"] = agg_df["Map_pred"].map(y_map)

    ax = plt.subplot(2, 3, 2)

    color_map = {"SIFT": "tab:blue", "PolyPhen": "tab:orange"}
    marker_map = {"SIFT": "o", "PolyPhen": "o"}

    for algo in agg_df["algorithm"].unique():
        sub = agg_df[agg_df["algorithm"] == algo]

        ax.scatter(
            sub["pos_bin"],
            sub["y"],
            s=sub["count"] * 8, 
            alpha=0.7,
            facecolors=color_map[algo],
            edgecolors="black",
            linewidths=0.3,
            marker=marker_map[algo],
            label=algo,
        )
    ax.set_yticks(range(len(impact_order)))
    ax.set_yticklabels(impact_order)
    ax.set_xlabel("Protein position (binned, 20 aa)")
    ax.set_title("B. Predicted functional effects across protein")
    ax.legend(title="Source", frameon=False, loc=0, bbox_to_anchor=(1, 1), borderaxespad=0.)
    ax.grid(True, alpha=0.3)
    # C
    plt.subplot(2, 3, 3)

    consequence_order = ['missense', 'frameshift', 'stop gained', '-', 'inframe deletion', 'insertion', 'stop lost']
    consequence_counts = df_variants['consequence'].value_counts()
    consequence_counts = consequence_counts.reindex(consequence_order, fill_value=0)
    counts = consequence_counts.copy()
    threshold = 0.03 * counts.sum()
    small = counts[counts < threshold].sum()
    counts = counts[counts >= threshold]
    counts['Other'] = small


    wedges, texts, autotexts = plt.pie(
        counts.values,
        labels=None,             
        autopct='%1.1f%%',        
        startangle=40,
        colors=plt.cm.tab20.colors
    )
    plt.legend(wedges, consequence_counts.index, title="Variant Consequences", bbox_to_anchor=(1.05, 0.5), loc="center left")
    plt.title('C. Distribution of Variant Consequences')
    plt.tight_layout()
    # D -- Plot PolyPhen scores only
    plt.subplot(2, 3, 4)

    # select PolyPhen predictions (case-insensitive)
    polyphen = polyphen_df[['position', 'score']].dropna()

    if polyphen.empty:
        plt.text(0.5, 0.5, 'No PolyPhen predictions available', ha='center', va='center')
        plt.xlabel('Amino acid position')
        plt.ylabel('PolyPhen score')
    else:
        vmin = polyphen['score'].min()
        vmax = polyphen['score'].max()
        sc = plt.scatter(
            polyphen['position'],
            polyphen['score'],
            c=polyphen['score'],
            cmap='GnBu',
            vmin=vmin,
            vmax=vmax,
            s=30,
            alpha=0.85
        )
        mappable = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=vmin, vmax=vmax), cmap='GnBu')
        mappable.set_array(polyphen['score'].values)
        ax = plt.gca()
        fig = plt.gcf()
        fig.colorbar(mappable, ax=ax, label=None)
        plt.xlabel('Amino acid position')
        plt.ylabel('PolyPhen score')
        plt.title('D. PolyPhen prediction scores across protein sequence')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
    # E
    plt.subplot(2, 3, 5)

    high_impact = polyphen_df[polyphen_df["impact_class"] == "High impact"]
    low_impact = polyphen_df[polyphen_df["impact_class"] == "Low / neutral"]

    plt.hist(
        [high_impact["position"], low_impact["position"]],
        color=["salmon", "darkblue"],
        bins=20,
        label=["High impact", "Low / neutral"],
        alpha=0.8
    )

    plt.xlabel("Amino acid position")
    plt.ylabel("Number of variants")
    plt.title("E. Positional enrichment of variant impact")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    # F
    plt.subplot(2, 3, 6)

    sns.countplot(
        data=pred_df,
        x="impact_class",
        hue="algorithm",
        order=["High impact", "Moderate impact", "Low / neutral", "Uncertain"],
        palette="Paired"
    )

    plt.xlabel("Impact class")
    plt.ylabel("Number of variants")
    plt.title("F. Predicted impact by algorithm")
    plt.xticks(rotation=30, ha="right")
    plt.legend(title="Predictor")
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    #plt.show()
    summary_text = f"""
    \n📈 VARIANT ANALYSIS SUMMARY
    Total variants: {df_variants['variant_id'].nunique()}
    Protein position range: {polyphen_df['position'].min()}–{polyphen_df['position'].max()}
    \nImpact class counts:
    {polyphen_df['impact_class'].value_counts()}
    \nPercent predicted high impact: {100 * len(polyphen_df[polyphen_df['impact_class'] == 'High impact']) / len(polyphen_df):.2f}%
    Top deleterious variants:
    {polyphen_df.sort_values(by='score', ascending=False).head(5)[['variant_id', 'position', 'score', 'prediction']].to_string(index=False)}
    \nVariants with disease association: {df_variants[df_variants['association'].notna()]['variant_id'].nunique()}

    """
    explain = """
    Plot A shows the distribution of variants along the protein. Variants cluster, with high-impact variants concentrated at positions critical for function.
    Plot B visualizes predicted functional effects by source (SIFT vs PolyPhen). Larger points indicate more variants in that bin.
    Plot C summarizes the types of variant consequences observed, with missense mutations being most common.
    Plot D displays PolyPhen prediction scores across the protein sequence, with higher scores indicating more likely damaging effects.
    Plot E highlights positional enrichment of high vs low impact variants, showing clustering of high-impact variants at key regions.
    Plot F compares predicted impact classes by algorithm, revealing differences in sensitivity between predictors.
    """
    return fig, summary_text, explain

#fig, summary = Variant_analysis("P61073")

def disease_associated_variants(uniprot_id):
    df_variants, pred_df = variant_dataframe(uniprot_id)
        # dISEASE ASSOCIATIONS EXTRACTION
    def extract_diseases(association):
        if not association:
            return None

        names = []
        if isinstance(association, list):
            for item in association:
                if isinstance(item, dict) and item.get('disease') is True:
                    if item.get('name'):
                        names.append(item['name'])

        return "; ".join(names) if names else None

    def extract_polyphen(predictions):
        if not predictions or not isinstance(predictions, list):
            return None, None

        for pred in predictions:
            if pred.get('predAlgorithmNameType') == 'PolyPhen':
                score = pred.get('score')
                label = pred.get('predictionValType')
                return score, label

        return None, None

    df_variants[['PolyPhen_score', 'PolyPhen_prediction']] = (
        df_variants['predictions']
        .apply(lambda x: pd.Series(extract_polyphen(x)))
    )
    df_variants['DiseaseList'] = df_variants['association'].apply(extract_diseases)

    df_disease_exploded = df_variants.explode('DiseaseList')
    df_disease_exploded['DiseaseList'] = df_disease_exploded['DiseaseList'].fillna('Not specified')

    disease_counts = df_disease_exploded['DiseaseList'].value_counts()

    total_variants = len(df_variants)

    summary = "Disease Distribution:\n"
    for disease, count in disease_counts.items():
        summary += f" • {disease}: {count} ({count/total_variants*100:.1f}%)\n"

    df_variants["begin"] = pd.to_numeric(df_variants["begin"], errors="coerce")
    protein_length = df_variants['begin'].max() + 10
    
    ### disease-associated variants table
    df_disease_table = df_variants.copy()
    df_disease_table['Disease'] = df_disease_table['association'].apply(extract_diseases)

    df_disease_table = df_disease_table[df_disease_table['Disease'].notnull()]

    df_disease_table = df_disease_table[[
        'variant_id',
        'external_url',
        'begin',
        'genomicLocation',
        'consequence',
        'wild_type',
        'mutatedType',
        'Disease',
        'PolyPhen_prediction'
    ]]

    df_disease_table = df_disease_table.sort_values('begin').reset_index(drop=True)

    ### disease plotting
    df_hotspot = df_variants[
    (df_variants["DiseaseList"].notna()) &
    (df_variants["begin"].notna())].copy()

    bin_size = 20
    bins = list(range(0, int(protein_length) + bin_size, bin_size))

    df_hotspot["position_bin"] = pd.cut(
        df_hotspot["begin"],
        bins=bins,
        include_lowest=True
    )

    hotspot_counts = (
        df_hotspot
        .groupby("position_bin")
        .size()
        .reset_index(name="DiseaseVariantCount")
    )

    hotspot_counts["bin_center"] = hotspot_counts["position_bin"].apply(
        lambda x: int((x.left + x.right) / 2)
    )

    df_variants["IsDamaging"] = (
        (df_variants["PolyPhen_prediction"] == "probably damaging") &
        (df_variants["PolyPhen_score"] >= 0.85)
    )

    df_priority = df_variants[
        (df_variants["IsDamaging"]) &
        (df_variants["DiseaseList"].notna())
    ].copy()

    total = len(df_variants)
    damaging = df_variants["IsDamaging"].sum()
    disease = df_variants["DiseaseList"].notna().sum()
    high_risk = len(df_priority)

    top_variants = (
        df_priority
        .sort_values("PolyPhen_score", ascending=False)
        .head(10)
    )

    top_variants[[
        "begin",
        "wild_type",
        "alt_seq",
        "PolyPhen_score",
        "DiseaseList"
    ]]

    df_diseases = df_priority.copy()

    disease_counts = (
        df_diseases["DiseaseList"]
        .str.split("; ")
        .explode()
        .value_counts()
        .head(10)
    )
    

    summary = pd.DataFrame({
        "Group": ["All variants", "Disease-associated", "Damaging", "High-risk"],
        "Count": [
            len(df_variants),
            disease,
            damaging,
            high_risk
        ]
    })

    text2 = f"""
    Disease Causing Variant Summary

    • Total variants: {total}

    • Predicted damaging variants:
    └─ {damaging} ({damaging/total*100:.1f}%)

    • Disease-associated variants:
    └─ {disease} ({disease/total*100:.1f}%)

    • High-confidence disease-damaging variants:
    └─ {high_risk} ({high_risk/total*100:.1f}%)

    Most frequent diseases linked to damaging variants:
    {disease_counts}
    """.strip()

    
    #--------Plotting---------
    sns.set_style("whitegrid")

    fig2 = plt.figure(figsize=(16, 10))

    gs = fig2.add_gridspec(
        2, 2,
        height_ratios=[1.1, 1],
        width_ratios=[0.7, 1.3],  # more space for C
        hspace=0.5,
        wspace=0.9
    )

    sns.set_style("whitegrid")

    #fig = plt.figure(figsize=(16, 10))

    # ------A. Disease variant distribution along sequence ----------------
    ax1 = fig2.add_subplot(gs[0, :])

    sns.barplot(
        data=hotspot_counts,
        x="bin_center",
        y="DiseaseVariantCount",
        color="#6baed6",
        edgecolor="black",
        ax=ax1
    )

    ax1.set_title("A. Distribution of Disease-Associated Variants Along Protein Sequence", fontsize=12)
    ax1.set_xlabel("Protein Position (amino acid, bin centers)")
    ax1.set_ylabel("Number of Disease-Associated Variants")
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3)

    # ---------------- B. Variant filtering summary  ----------------
    ax2 = fig2.add_subplot(gs[1,0])

    colors_b = sns.color_palette("PuBu", n_colors=len(summary))


    wedges, texts, autotexts = ax2.pie(
        summary["Count"],
        labels=None,                      
        autopct=lambda p: f"{p:.1f}%" if p > 4 else "",
        startangle=90,
        colors=colors_b,
        pctdistance=0.7,
        wedgeprops=dict(edgecolor="white")
    )

    ax2.legend(
        wedges,
        summary["Group"],
        title="Variant category",
        bbox_to_anchor=(1.05, 0.5), 
        loc="upper center",
        frameon=False,
        fontsize=6,
        title_fontsize=8
    )
    for autotext in autotexts:
        autotext.set_fontsize(7)

    ax2.set_title("B. Variant Filtering and Prioritization Summary",fontsize=11,pad=10)
    ax2.axis("equal")
    # ---------------- C. Disease frequency ----------------
    ax3 = fig2.add_subplot(gs[1,1])

    wrapped_labels = [
        "\n".join(textwrap.wrap(d, 35))
        for d in disease_counts.index
    ]

    sns.barplot(
        x=disease_counts.values,
        y=wrapped_labels,
        palette=sns.color_palette("Blues_d", len(wrapped_labels)),
        edgecolor="black",
        ax=ax3
    )

    ax3.yaxis.set_label_coords(-0.18, 0.5)
    ax3.tick_params(axis='y', labelsize=7)
    ax3.tick_params(axis='x', labelsize=8)
    ax3.set_title("C. Most Frequent Diseases Associated with\nDamaging Variants",fontsize=10)
    ax3.set_xlabel("Number of High-Risk Variants")
    #ax3.yaxis.set_label_coords(-0.12, 0.5)

    plt.tight_layout()

    #print("\nDisease-associated Variants Table:")
    return df_disease_table, summary, text2, fig2


#disease_associated_variants('Q9Y243')
