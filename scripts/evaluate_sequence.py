# scripts/evaluate_sequences.py
# Threading de secuencias MPNN sobre backbone WT 1DKL + FastRelax (ref2015)
# + scoring total. Input: top10_queries.fasta y fitasa_backbone_clean_1DKL.pdb

import os
from pyrosetta import init, pose_from_pdb, Pose, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

# Inicializar PyRosetta con rotámeros extendidos
init("-ex1 -ex2 -use_input_sc")

# Rutas (ajustar al entorno del lab)
backbone_path = "./input/fitasa_backbone_clean_1DKL.pdb"
fasta_path    = "./input/top10_queries.fasta"
out_dir       = "./output/threaded_structures"
scores_csv    = "./output/scores.csv"

pose_template = pose_from_pdb(backbone_path)

with open(fasta_path, 'r') as f:
    lines = f.read().splitlines()
    seq_headers = lines[::2]
    sequences = lines[1::2]

os.makedirs(out_dir, exist_ok=True)

scorefxn = create_score_function("ref2015")
relax = FastRelax()
relax.set_scorefxn(scorefxn)


def thread_sequence(base_pose: Pose, sequence: str, header: str) -> Pose:
    sequence = sequence.strip()
    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")

    if not set(sequence).issubset(valid_aas):
        raise ValueError(f"Secuencia inválida en {header}: contiene caracteres fuera de los 20 AAs estándar.")

    if len(sequence) != base_pose.total_residue():
        raise ValueError(f"Longitud de secuencia y PDB no coinciden en {header}: {len(sequence)} vs {base_pose.total_residue()}")

    print(f"Secuencia válida {header} ({len(sequence)} aa)")

    threaded = Pose()
    threaded.assign(base_pose)
    for i, aa in enumerate(sequence, start=1):
        if aa != threaded.residue(i).name1():
            MutateResidue(i, aa).apply(threaded)
    return threaded


assert len(seq_headers) == len(sequences), "FASTA mal formateado"

results = []
for header, sequence in zip(seq_headers, sequences):
    print(f"Procesando {header}...")
    pose = thread_sequence(pose_template, sequence.strip(), header)
    relax.apply(pose)
    total_score = scorefxn(pose)

    out_pdb = f"{out_dir}/{header[1:]}.pdb"
    pose.dump_pdb(out_pdb)
    results.append((header[1:], total_score))

with open(scores_csv, 'w') as f:
    f.write("ID,Total_REU\n")
    for name, score in sorted(results, key=lambda x: x[1]):
        f.write(f"{name},{score:.3f}\n")

print("Evaluación completa.")
