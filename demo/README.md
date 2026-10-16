# Small, complete synthetic examples

All inputs are included. These are deterministic synthetic installation examples, **not biological validation data**.

After installing and activating the environment:

```bash
bash demo/run_all.sh --cores 2
```

| Strategy | Input | Required result |
| --- | --- | --- |
| I | 50-kb reference and query containing an inserted 360-bp sequence | Recover the 360-bp InDel and produce at least one usable primer pair |
| II | 12-kb reference, paired 150-bp reads from a reference parent and a parent with a 40-bp deletion | Recover the deletion with `p1=0/0`, `p2=1/1`, and produce primers |
| III | 12-kb genomes, one 1.2-kb CDS, and a 30-bp insertion in the second genome | Produce an alignment containing the CDS and both genomes, with the known size difference |

The paired reads tile variable-length fragments with distinct coordinates and matching R1/R2 identifiers, providing enough coverage to exercise duplicate marking and variant calling.

`generate_demo_inputs.py` rebuilds all inputs reproducibly with fixed random seeds and only the Python standard library. The runner calls it automatically and leaves unchanged files untouched. `expected.json` records the known variants and genotype expectations.

The assembly example includes a downstream SNP as well as the target insertion. SyRI 1.8.2 can omit an isolated final InDel when combining its output; the known-InDel check guards against an empty result. The example is an installation check, not an exhaustive assessment of the variant caller's sensitivity.

Results are written under `demo/work/strategy1/`, `strategy2/`, and `strategy3/`. Read `validation.json` in each output directory for the content-check result. Strategy I and II primer designs are in `primer3/primer_markers.xlsx`; Strategy III alignments are listed in `strategy3/alignment/manifest.tsv`.

Use `bash demo/run_all.sh --strategy 2 --cores 2` to run only the read-based example. No input downloads are required after cloning or extracting the ZIP.
