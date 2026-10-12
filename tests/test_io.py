from dataclasses import replace
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from motif_audit import cli
from motif_audit.core import AuditError, Budget, analyze, digest
from motif_audit.report import match_tsv, render
from test_core import manifest

ROOT = Path(__file__).resolve().parents[1]


class OutputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.output = Path(self.temp.name) / "output"
        self.result = analyze(b">=hostile </script><img src=https://invalid.test onerror=alert(1)> @@JS@@\nACGTACATAT", manifest("ACGTAC", "ATAT"), shuffles=2)

    def test_complete_artifact_hashes_and_inert_html(self):
        cli.write_output(self.output, self.result, Budget())
        complete = json.loads((self.output / "completion.json").read_text())
        self.assertEqual(set(p.name for p in self.output.iterdir()), {"completion.json", "analysis.json", "matches.tsv", "report.html"})
        for name, meta in complete["artifacts"].items():
            content = (self.output / name).read_bytes()
            self.assertEqual(len(content), meta["bytes"])
            self.assertEqual(digest(content), meta["sha256"])
        document = (self.output / "report.html").read_text()
        self.assertNotIn("<img src=", document)
        self.assertIn("\\u003c/script\\u003e", document)
        self.assertEqual(document.count("<script"), 2)
        self.assertIn("connect-src 'none'", document)
        self.assertIn("'"+"=hostile", (self.output / "matches.tsv").read_text())

    def test_collision_keeps_existing_data(self):
        self.output.mkdir()
        sentinel = self.output / "sentinel"
        sentinel.write_text("keep")
        with self.assertRaises(FileExistsError):
            cli.write_output(self.output, self.result, Budget())
        self.assertEqual(sentinel.read_text(), "keep")
        sentinel.unlink(); self.output.rmdir()
        self.output.symlink_to(Path(self.temp.name), target_is_directory=True)
        with self.assertRaises(FileExistsError):
            cli.write_output(self.output, self.result, Budget())
        self.assertTrue(self.output.is_symlink())

    def test_serialization_size_cancellation_and_write_cleanup(self):
        for function, error in [("canonical", TypeError("bad serialization")),
                                ("render", ValueError("bad template")),
                                ("match_tsv", KeyboardInterrupt())]:
            with self.subTest(function=function):
                with patch.object(cli, function, side_effect=error):
                    with self.assertRaises(type(error)):
                        cli.write_output(self.output, self.result, Budget())
                self.assertFalse(self.output.exists())
        with patch("motif_audit.cli.os.fsync", side_effect=OSError("disk full")):
            with self.assertRaises(OSError):
                cli.write_output(self.output, self.result, Budget())
        self.assertFalse(self.output.exists())
        with patch("pathlib.Path.rename", side_effect=OSError("rename failure")):
            with self.assertRaises(OSError):
                cli.write_output(self.output, self.result, Budget())
        self.assertFalse(self.output.exists())
        with self.assertRaises(AuditError):
            cli.write_output(self.output, self.result, Budget(-1))
        self.assertFalse(self.output.exists())
        self.result["settings"]["limits"]["artifact_bytes"] = 10
        with self.assertRaises(AuditError):
            cli.write_output(self.output, self.result, Budget())
        self.assertFalse(self.output.exists())

    def test_empty_results_export(self):
        result = analyze(b">x\nNNNN", manifest("A"), shuffles=2)
        self.assertEqual(len(match_tsv(result).splitlines()), 1)
        self.assertIn(b"<html", render(result))
        self.assertEqual(result["statistics"][0]["observed"], 0)
        self.assertEqual(result["statistics"][0]["p_holm"], 1)

    def test_special_file_and_read_byte_limit(self):
        fifo = Path(self.temp.name) / "fifo"
        os.mkfifo(fifo)
        with self.assertRaises(AuditError):
            cli.read_bounded(fifo, 10)
        f = Path(self.temp.name) / "big"
        f.write_bytes(b"A"*11)
        with self.assertRaises(AuditError):
            cli.read_bounded(f, 10)

    def test_cli_validation_and_deadline(self):
        base = [sys.executable, "-m", "motif_audit", "--fasta", str(ROOT / "examples/sequences.fasta"),
                "--motifs", str(ROOT / "examples/motifs.json"), "--out", str(self.output)]
        for args in [["--seconds", "nan"], ["--seconds", "inf"], ["--seconds", "0"],
                     ["--seconds", "121"], ["--seconds", "0.000001"], ["--shuffles", "0"], ["--seed", "-1"]]:
            p = subprocess.run(base + args, capture_output=True, cwd=ROOT, timeout=10)
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertFalse(self.output.exists())

    def test_actual_signals_cleanup(self):
        # Synchronize at a known point AFTER output reservation, before serialization.
        # This avoids a flaky timing race between Python startup and signal handlers.
        for sig in [signal.SIGINT, signal.SIGTERM, signal.SIGALRM]:
            code = '''
from pathlib import Path
import sys, time
from motif_audit import cli
from motif_audit.core import Budget
out = Path(sys.argv[1])
def slow(result):
    (out / 'ready').write_text('ready')
    time.sleep(10)
cli.canonical = slow
try:
    with cli.deadline(20):
        cli.write_output(out, {}, Budget())
except BaseException:
    sys.exit(9)
'''
            process = subprocess.Popen([sys.executable, "-c", code, str(self.output)], cwd=ROOT)
            try:
                deadline = time.monotonic()+5
                while not (self.output / "ready").exists() and time.monotonic() < deadline:
                    time.sleep(.01)
                self.assertTrue((self.output / "ready").exists())
                process.send_signal(sig)
                self.assertEqual(process.wait(timeout=5), 9)
                self.assertFalse(self.output.exists())
            finally:
                if process.poll() is None:
                    process.kill(); process.wait()


if __name__ == "__main__":
    unittest.main()
