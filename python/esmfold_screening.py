import os
import argparse
import torch
import string
from Bio.SeqIO import FastaIO
import esm
from esm.esmfold.v1.esmfold import ESMFold

def clean_sequence(seq):
    """Removes non-standard amino acids that crash ESMFold"""
    allowed = set(string.ascii_uppercase) - set("BJOUXZ")
    return "".join([c for c in seq.upper() if c in allowed])


def main():
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description="Run 1-vs-All ESMFold cofolding screening on Rorqual.")
    parser.add_argument("--bait", type=str, required=True, help="Amino acid sequence of the bait protein.")
    parser.add_argument("--fasta", type=str, required=True, help="Path to the Arabidopsis proteome multi-FASTA file.")
    parser.add_argument("--outdir", type=str, default="screening_results",
                        help="Directory to save generated PDB files.")
    parser.add_argument("--chunk-size", type=int, default=None,
                        help="Axial attention chunk size for ESMFold (model.set_chunk_size). "
                             "Lower values reduce peak GPU memory at the cost of speed. "
                             "Omit for no chunking (short sequences only).")

    args = parser.parse_args()

    # Process inputs
    bait_seq = clean_sequence(args.bait)
    fasta_path = args.fasta
    output_dir = args.outdir

    os.makedirs(output_dir, exist_ok=True)

    print("Loading ESMFold on H100 GPU...")
    original_esm2_loader = esm.pretrained.esm2_t36_3B_UR50D
    esm.pretrained.esm2_t36_3B_UR50D = lambda: esm.esmfold.v1.misc.load_esm_model(
        esm_type="esm2_t36_3B_UR50D",
        preload_weights=False  # <--- CRITICAL: Stops the internet download trigger
    )

    weight_path = "/home/jflucier/links/scratch/programs/esm/esmfold_3B_v1.pt"
    model_data = torch.load(weight_path, map_location="cpu", weights_only=False)

    # 2. Extract the model architecture configurations embedded inside the file
    cfg = model_data["cfg"]["model"]

    # 3. Instantiate the structural network matrix completely locally
    model = ESMFold(esmfold_config=cfg)

    # 4. Bind the offline weights to the new architecture shell
    model.load_state_dict(model_data["model"])

    # 5. Push the compiled model cleanly onto your H100 MIG slice
    model = model.eval().cuda()

    # Optimize matrix operations specifically for the H100 SXM5 architecture
    torch.set_float32_matmul_precision('high')

    if args.chunk_size is not None:
        print(f"Enabling axial attention chunking with chunk_size={args.chunk_size}")
        model.set_chunk_size(args.chunk_size)

    print(f"Starting screen using bait length: {len(bait_seq)} AA")
    print(f"Reading targets from: {fasta_path}")
    print(f"Saving structures to: {output_dir}")

    # Stream through the large proteome file to conserve RAM
    with open(fasta_path, "r") as handle:
        for record in FastaIO.SimpleFastaParser(handle):
            header_line = record[0]
            target_seq = clean_sequence(record[1])

            # Extract first word of the header as the unique filename token (e.g., AT1G01010)
            target_id = header_line.split()[0].replace("|", "_").replace("/", "_")

            # Filter out empty entries or excessive sequence lengths (>1200 total residues can trigger OOMs)
            # if not target_seq or (len(bait_seq) + len(target_seq) > 1200):
            #     continue

            output_pdb = os.path.join(output_dir, f"{target_id}.pdb")
            if os.path.exists(output_pdb):
                continue  # Skips execution if this pairing was already finished

            # Colon ':' triggers multi-chain cofolding inside ESMFold
            cofold_sequence = f"{bait_seq}:{target_seq}"

            try:
                with torch.no_grad():
                    output_string = model.infer_pdb(cofold_sequence)

                    with open(output_pdb, "w") as f:
                        f.write(output_string)

                print(f"Successfully cofolded: bait vs {target_id}")

            except Exception as e:
                print(f"Skipping {target_id} due to an execution error: {e}")
                # Clear cached (but unused) CUDA memory so a single OOM/error doesn't
                # cause subsequent, otherwise-fittable targets to also fail.
                torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
