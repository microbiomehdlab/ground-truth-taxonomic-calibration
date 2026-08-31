import gzip
import pathlib
import subprocess
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class PreproductionIntegrityTests(unittest.TestCase):
    def test_pure_pool_audit_is_isolated_and_deterministic(self):
        runner = (ROOT / "datasets/yachida/run_pure_pool_audit.sh").read_text()
        sbatch = (ROOT / "run_yachida_pure_pool_audit.sbatch").read_text()
        self.assertIn('PURE_POOL_AUDIT_PAIRS="${PURE_POOL_AUDIT_PAIRS:-1000000}"', runner)
        self.assertIn("selection_rule\\tdeterministic_prefix", runner)
        self.assertIn('PURE_POOL_AUDIT_ROOT', runner)
        self.assertNotIn('PERSISTENT_RESULTS_ROOT', runner)
        self.assertIn("#SBATCH --cpus-per-task=8", sbatch)
        self.assertIn("#SBATCH --mem=48G", sbatch)

    def test_yachida_production_resources_fit_standard_nodes(self):
        script = (ROOT / "run_yachida_batch_sample.sbatch").read_text()
        self.assertIn("#SBATCH --cpus-per-task=16", script)
        self.assertIn("#SBATCH --mem=48G", script)
        runner = (ROOT / "datasets/yachida/run_full_streaming_sample.sh").read_text()
        self.assertIn('PROFILE_CONCURRENCY="${PROFILE_CONCURRENCY:-1}"', runner)

    def test_preproduction_gate_clears_stale_success_marker(self):
        script = (ROOT / "run_yachida_preproduction_audit.sh").read_text()
        marker_cleanup = 'rm -f -- "$AUDIT_DIR/SUCCESS"'
        self.assertIn(marker_cleanup, script)
        self.assertLess(script.index(marker_cleanup), script.index("python3 \"$ROOT/datasets/yachida/preproduction_audit.py\""))

    def test_generic_bracken_threshold_is_not_cpu_thread_count(self):
        script = (ROOT / "workflows/kraken2_bracken/classify_bracken.sbatch").read_text()
        self.assertIn('BRACKEN_THRESHOLD="${BRACKEN_THRESHOLD:-10}"', script)
        self.assertIn('-t "$BRACKEN_THRESHOLD"', script)
        self.assertNotIn('-t "$THREADS"', script)

    def test_community_allocation_is_exact_and_nonnegative(self):
        script = ROOT / "spikes/scripts/spikein/allocate_community_reads.py"
        for total in range(1, 101):
            command = ["python3", str(script), "--total", str(total)]
            for label in "ABCDEFGHIJ":
                command.extend(("--member", f"{label}:1"))
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            counts = [int(line.split("\t")[1]) for line in result.stdout.splitlines()[1:]]
            self.assertEqual(sum(counts), total)
            self.assertGreaterEqual(min(counts), 0)

    def test_paired_fastq_validator_accepts_synchronized_mates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            r1, r2, output = root / "r1.fq.gz", root / "r2.fq.gz", root / "integrity.tsv"
            with gzip.open(r1, "wt") as left, gzip.open(r2, "wt") as right:
                left.write("@readA/1\nACGT\n+\nIIII\n@readB 1:N:0:1\nGGCC\n+\nIIII\n")
                right.write("@readA/2\nTGCA\n+\nIIII\n@readB 2:N:0:1\nCCGG\n+\nIIII\n")
            subprocess.run(
                ["python3", str(ROOT / "scripts/validate_paired_fastq.py"), "--r1", str(r1), "--r2", str(r2), "--output", str(output)],
                check=True,
            )
            self.assertIn("synchronized_pairs\t2", output.read_text())

    def test_paired_fastq_validator_rejects_mismatched_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            r1, r2, output = root / "r1.fq", root / "r2.fq", root / "integrity.tsv"
            r1.write_text("@readA/1\nACGT\n+\nIIII\n")
            r2.write_text("@readB/2\nTGCA\n+\nIIII\n")
            result = subprocess.run(
                ["python3", str(ROOT / "scripts/validate_paired_fastq.py"), "--r1", str(r1), "--r2", str(r2), "--output", str(output)],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("mate identifiers differ", result.stderr)


if __name__ == "__main__":
    unittest.main()
