from Bio.PDB import MMCIFParser, PDBParser, PPBuilder, Superimposer
from Bio.PDB.vectors import calc_dihedral
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
import requests
import io
import numpy as np
import matplotlib.pyplot as plt
from requests.exceptions import HTTPError, RequestException

# ---------- INPUT ----------
#pdb_id     = "4PED"        # experimental structure
#uniprot_id = "P42336" #"Q96D53"      # UniProt ID
#chain_id   = "A"
##start_res, end_res = 1, 544
high_cut = 2.0

ppb = PPBuilder()

def pick_longest_structure(uniprot_id: str):
    
    uniprot_id = uniprot_id.strip().upper()
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    r = requests.get(url, timeout=30)
    
    # PDBe uses 404 to mean "no data for this UniProt"
    if r.status_code == 404:
        return None, None, None, None, None
    
    r.raise_for_status()
    data = r.json()

    candidates = data.get(uniprot_id, [])
    if not candidates:
        return None, None, None, None, None

    candidates = data.get(uniprot_id, [])
    if not candidates:
        return None, None, None, None, None  # no experimental structure found

    # prefer longest UniProt span; tie-break by best (lowest) resolution if available
    def span_len(x):
        return int(x.get("unp_end", 0)) - int(x.get("unp_start", 0)) + 1

    def res_value(x):
        return x["resolution"] if x.get("resolution") is not None else 99.0

    best = sorted(
        candidates,
        key=lambda x: (span_len(x), -1.0 / res_value(x)),  # longer then lower resolution
        reverse=True
    )[0]

    pdb_id = best["pdb_id"].upper()
    chain_id = best["chain_id"]
    start_res = int(best["unp_start"])
    end_res = int(best["unp_end"])
    return pdb_id, chain_id, start_res, end_res, best



# 1. FETCH REMOTELY

def fetch_pdb_mmcif(pdb_id):
    """Fetch experimental structure as mmCIF from RCSB"""
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.cif"
    print(f"[INFO] Fetching PDB {pdb_id} from: {url}")
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        return io.StringIO(r.text)
    except HTTPError as e:
        print(f"[ERROR] Failed to fetch PDB {pdb_id}: {e}")
        raise

def fetch_alphafold_pdb(uniprot_id):
    """Fetch AlphaFold PDB using the API to find the latest file URL"""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    print(f"[INFO] Querying AlphaFold API for: {uniprot_id}")
    try:
        r = requests.get(api_url, timeout=30)
        r.raise_for_status()
        data = r.json()
        if not data:
            raise ValueError(f"No AlphaFold entry found for {uniprot_id}")
        
        pdb_url = data[0]['pdbUrl']
        print(f"[INFO] Downloading PDB from: {pdb_url}")
        
        pdb_r = requests.get(pdb_url, timeout=60)
        pdb_r.raise_for_status()
        return io.StringIO(pdb_r.text)
    except (HTTPError, KeyError, IndexError, ValueError) as e:
        print(f"[ERROR] AlphaFold fetch failed: {e}")
        raise

def fetch_uniprot_fasta(uniprot_id):
    """Fetch UniProt reference sequence"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        handle = io.StringIO(r.text)
        record = SeqIO.read(handle, "fasta")
        return str(record.seq)
    except HTTPError as e:
        print(f"[ERROR] UniProt {uniprot_id} not found: {e}")
        raise

# Structure Comparison Helper Functions
def get_ca_range(chain, start, end):
    ca_atoms = []
    for res in chain:
        n = res.id[1]
        if start <= n <= end and "CA" in res:
            ca_atoms.append(res["CA"])
    return ca_atoms

def per_residue_rmsd_range(your_chain, af_chain, start, end):
    your_res = {r.id[1]: r for r in your_chain.get_residues() if "CA" in r}
    af_res   = {r.id[1]: r for r in af_chain.get_residues() if "CA" in r}
    
    common_nums = sorted(set(your_res.keys()) & set(af_res.keys()) & set(range(start, end+1)))
    
    rmsd_vals = []
    res_nums  = []
    for n in common_nums:
        diff = your_res[n]["CA"] - af_res[n]["CA"]
        rmsd_vals.append(np.sqrt(diff * diff))
        res_nums.append(n)
    return res_nums, rmsd_vals

def segment_rmsd(your_chain, af_chain, seg_dict):
    seg_stats = {}
    for name, (s, e) in seg_dict.items():
        res_nums_seg, rmsd_vals = per_residue_rmsd_range(your_chain, af_chain, s, e)
        if len(rmsd_vals) == 0: continue
        seg_stats[name] = {
            "mean": float(np.mean(rmsd_vals)),
            "std":  float(np.std(rmsd_vals)),
            "max":  float(np.max(rmsd_vals)),
        }
    return seg_stats

def ca_distance_matrix(ca_atoms):
    n = len(ca_atoms)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = ca_atoms[i] - ca_atoms[j]
            mat[i, j] = mat[j, i] = np.sqrt(d*d)
    return mat


# 2. PROTEIN IDENTITY VERIFICATION
def get_full_sequence(chain):
    """Joins all peptide fragments to get the full sequence present in the structure"""
    peptides = ppb.build_peptides(chain)
    return "".join([str(pp.get_sequence()) for pp in peptides])

def verify_protein_identity(struct_chain, ref_seq):
    """Calculates sequence identity between the structure and reference sequence"""
    struct_seq = get_full_sequence(struct_chain)
    
    # Perform alignment
    matrix = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globalds(struct_seq, ref_seq, matrix, -10, -0.5)
    
    if not alignments:
        return {
            "struct_length": len(struct_seq),
            "ref_length": len(ref_seq),
            "identity": 0.0,
            "is_same_protein": False
        }
    
    best_aln = alignments[0]
    matches = sum(1 for a, b in zip(best_aln.seqA, best_aln.seqB) if a == b)
    identity = (matches / len(ref_seq)) * 100
    
    return {
        "struct_length": len(struct_seq),
        "ref_length": len(ref_seq),
        "identity": identity,
        "is_same_protein": identity > 90.0
    }

def run_full_verification(pdb_id, uniprot_id, chain_id="A"):
    # Fetch structures
    cif_handle = fetch_pdb_mmcif(pdb_id)
    af_handle  = fetch_alphafold_pdb(uniprot_id)

    # Parse structures
    cif_parser = MMCIFParser(QUIET=True)
    pdb_parser = PDBParser(QUIET=True)

    your_struct = cif_parser.get_structure("your", cif_handle)
    af_struct   = pdb_parser.get_structure("af",   af_handle)

    # Correctly identify chains
    your_model = list(your_struct.get_models())[0]
    your_chain = your_model[chain_id]
    
    af_model = list(af_struct.get_models())[0]
    # AF files usually only have one chain, often named 'A'
    if chain_id in af_model:
        af_chain = af_model[chain_id]
    else:
        af_chain = list(af_model.get_chains())[0]

    # Verify against UniProt
    ref_seq = fetch_uniprot_fasta(uniprot_id)
    your_result = verify_protein_identity(your_chain, ref_seq)
    af_result   = verify_protein_identity(af_chain,   ref_seq)

    summary = f"""
    PROTEIN IDENTITY VERIFICATION
    {'='*70}
    PDB ID:                 {pdb_id}
    UniProt ID:             {uniprot_id}
    Experimental length:    {your_result['struct_length']} aa
    AlphaFold length:       {af_result['struct_length']} aa
    UniProt reference:      {your_result['ref_length']} aa
    Exp vs UniProt identity:{your_result['identity']:.2f}%
    AF vs UniProt identity: {af_result['identity']:.2f}%

    FINAL VERDICT
    {'-'*70}
    Experimental: {'✅ MATCH' if your_result['is_same_protein'] else '❌ MISMATCH'}
    AlphaFold:    {'✅ MATCH' if af_result['is_same_protein'] else '❌ MISMATCH'}
    {'='*70}
    """

    return summary,your_result, af_result, your_chain, af_chain



# ---
def structural_comparison(uniprot_id):
    
    pdb_id, chain_id, start_res, end_res, meta = pick_longest_structure(uniprot_id)

    if pdb_id is None:
        fail_summary = f"""
        No experimental PDB found for {uniprot_id}. Consider AlphaFold fallback.
        PROTEIN IDENTITY VERIFICATION
        {'='*70}
        UniProt ID:             {uniprot_id}
        Experimental structure: NOT AVAILABLE
                FINAL VERDICT
                {'-'*70}
                Experimental: ❌ NOT AVAILABLE
                
            """
        return None, fail_summary, None
    
    # Run verification
    summary, your_res_info, af_res_info, your_chain, af_chain = run_full_verification(pdb_id, uniprot_id, chain_id)

    # 1. Structural Superimposition
    your_ca = get_ca_range(your_chain, start_res, end_res)
    af_ca   = get_ca_range(af_chain,   start_res, end_res)
    
    your_ca_dict = {a.get_parent().id[1]: a for a in your_ca}
    af_ca_dict   = {a.get_parent().id[1]: a for a in af_ca}
    common_ids = sorted(set(your_ca_dict.keys()) & set(af_ca_dict.keys()))
    
    fixed_ca = [af_ca_dict[i] for i in common_ids]
    moving_ca = [your_ca_dict[i] for i in common_ids]

    sup = Superimposer()
    sup.set_atoms(fixed_ca, moving_ca)
    sup.apply(your_chain.get_atoms())

    # 2. Analysis
    res_nums, rmsd_per_res = per_residue_rmsd_range(your_chain, af_chain, start_res, end_res)
    
    your_ca_common = [your_ca_dict[i] for i in common_ids]
    af_ca_common   = [af_ca_dict[i] for i in common_ids]
    your_mat = ca_distance_matrix(your_ca_common)
    af_mat   = ca_distance_matrix(af_ca_common)

    # 3. Plotting
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # RMSD Plot
    axes[0,0].plot(res_nums, rmsd_per_res, color='red', lw=1)
    axes[0,0].axhline(high_cut, ls='--', color='orange')
    axes[0,0].set_title("Per-Residue RMSD")
    axes[0,0].set_ylabel("Å")

    # RMSD Histogram
    axes[0,1].hist(rmsd_per_res, bins=25, color='skyblue', edgecolor='black')
    axes[0,1].set_title("RMSD Distribution")

    # Distance Matrix Difference
    im = axes[1,0].imshow(np.abs(your_mat - af_mat), cmap='viridis')
    axes[1,0].set_title("Distance Matrix Difference")
    plt.colorbar(im, ax=axes[1,0])

    # Segment Bar Chart
    segments = {"N-Term": (1, 150), "Core": (151, 400), "C-Term": (401, 544)}
    stats = segment_rmsd(your_chain, af_chain, segments)
    axes[1,1].bar(stats.keys(), [s['mean'] for s in stats.values()], color='green')
    axes[1,1].set_title("Mean RMSD by Segment")

    plt.tight_layout()

    text = f"Global RMSD (on {len(common_ids)} Cα atoms): {sup.rms:.3f} Å"
    
    return fig, summary, text

#fig, summary, text = structural_comparison("Q96D53")
#print(summary)

