import os
import shutil
from pathlib import Path

# Source directories
TRENCH_DIR = Path(r"C:\Users\hooki\Trench")
VAULT_DIR = Path(r"C:\Users\hooki\trench_re_vault")

# Destination
CORPUS_DIR = TRENCH_DIR / "docs" / "notebooklm" / "staging_corpus"

def collect_files():
    if CORPUS_DIR.exists():
        shutil.rmtree(CORPUS_DIR)
        
    # Set up the internal folder structure
    dirs = {
        "re_notes": CORPUS_DIR / "reverse_engineering_notes",
        "official_theory": CORPUS_DIR / "theory_and_tutorials",
        "architecture": CORPUS_DIR / "architecture_and_math",
        "design_specs": CORPUS_DIR / "design_specs",
    }
    
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
        
    files_to_copy = [
        # RE Notes (from Vault)
        (VAULT_DIR / "analysis" / "morph_designer_decompile.md", dirs["re_notes"]),
        (VAULT_DIR / "analysis" / "contradictions.md", dirs["re_notes"]),
        (VAULT_DIR / "analysis" / "test_targets.md", dirs["re_notes"]),
        (VAULT_DIR / "analysis" / "open_questions.md", dirs["re_notes"]),
        (VAULT_DIR / "contracts" / "CPhantomMorphLPX" / "2026-03-16_MorphLP_zero_table.md", dirs["re_notes"]),
        (VAULT_DIR / "contracts" / "CPhantomMorphLPX" / "2026-03-16_comparison_MorphLP_vs_LPX.md", dirs["re_notes"]),
        (VAULT_DIR / "contracts" / "CPhantomMorphLPX" / "2026-03-16_vtable_decompilation.md", dirs["re_notes"]),
        
        # RE Notes (from Trench)
        (TRENCH_DIR / "docs" / "looperator" / "phantom_morph_lpx_contracts.md", dirs["re_notes"]),
        (TRENCH_DIR / "docs" / "looperator" / "talking_hedz_x3_surfaces.md", dirs["re_notes"]),
        (TRENCH_DIR / "dev" / "sift" / "lattice_law_evidence.md", dirs["re_notes"]),
        (TRENCH_DIR / "dev" / "sift" / "bodies_evidence.md", dirs["re_notes"]),
        (TRENCH_DIR / "docs" / "notebooklm" / "packs" / "ground-truth" / "trench_re_codex.md", dirs["re_notes"]),
        (TRENCH_DIR / "docs" / "notebooklm" / "packs" / "ground-truth" / "re-findings.md", dirs["re_notes"]),

        # Architecture & Math
        (TRENCH_DIR / "docs" / "notebooklm" / "packs" / "ground-truth" / "morpheus_cube_architecture_proven.md", dirs["architecture"]),
        (TRENCH_DIR / "docs" / "notebooklm" / "packs" / "ground-truth" / "q_law_firmware_proven.md", dirs["architecture"]),
        (TRENCH_DIR / "docs" / "notebooklm" / "packs" / "ground-truth" / "runtime_pipeline_proven.md", dirs["architecture"]),
        (TRENCH_DIR / "docs" / "notebooklm" / "packs" / "ground-truth" / "zero_trajectory_proven.md", dirs["architecture"]),
        (TRENCH_DIR / "docs" / "dsp-spec.md", dirs["architecture"]),
        
        # Theory & Tutorials
        (TRENCH_DIR / "docs" / "heritage_dillusionman_peak_shelf_morph.md", dirs["official_theory"]),
        
        # Design Specs & Targets
        (TRENCH_DIR / "BODIES.md", dirs["design_specs"]),
        (TRENCH_DIR / "SPEC.md", dirs["design_specs"]),
        (TRENCH_DIR / "docs" / "hostile_authoring_workflow_spec.md", dirs["design_specs"]),
        (TRENCH_DIR / "docs" / "world_q_behavior_gate_spec.md", dirs["design_specs"]),
    ]
    
    extracted_count = 0
    for src, dst_folder in files_to_copy:
        if src.exists():
            shutil.copy2(src, dst_folder / src.name)
            extracted_count += 1
            print(f"Copied: {src.name} -> {dst_folder.name}/")
        else:
            print(f"Skipping (not found): {src}")
            
    # Create the cookbook & glossary templates the user requested
    templates_dir = CORPUS_DIR / "templates_to_fill"
    templates_dir.mkdir(exist_ok=True)
    
    glossary_content = "# Parameter Glossary\n\n- **Morph**: ...\n- **Freq**: ...\n- **Q/Peak**: ...\n- **Shelf**: ...\n"
    (templates_dir / "parameter_glossary_template.md").write_text(glossary_content)
    
    cookbook_content = "# Design Cookbook\n\n- **Sub cleanup**\n- **Vocal push**\n- **Reese sweep**\n- **Metallic percussion**\n"
    (templates_dir / "patch_cookbook_template.md").write_text(cookbook_content)
    
    print(f"\nSuccessfully extracted {extracted_count} RE and design documents to {CORPUS_DIR}")

if __name__ == "__main__":
    collect_files()
