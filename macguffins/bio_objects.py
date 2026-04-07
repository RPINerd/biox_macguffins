"""All the classes that are used across the biox_macguffins codebase"""


import subprocess
from pathlib import Path

from macguffins.utils import extract_sample_info


class Primer:

    """A class to represent a primer"""

    def __init__(self, chr: str | int, num: int, start: int, end: int) -> None:
        """Initializes a primer object"""
        self.chr = chr
        self.number = num
        self.start = start
        self.end = end

    def __str__(self) -> str:
        """Returns a string representation of the primer"""
        return f"Primer {self.number} on chromosome {self.chr} from {self.start} to {self.end}"

    def __repr__(self) -> str:
        """Returns a string representation of the primer"""
        return f"Primer{self.number}, {self.chr}:{self.start}-{self.end}"


class RunCollection:

    """A class to represent a set of fastqs that are intended to be passed through a pipeline together"""

    def __init__(self, dir: Path, name: str = "") -> None:
        """Initializes a RunCollection object"""
        self.run_dir = dir
        self.run_name = name if name else dir.name
        self.samples: dict[str, RunSet] = {}
        self.collect_samples(dir)

    def __str__(self) -> str:
        """"""
        return self.run_name

    def __repr__(self) -> str:
        """"""
        return f"RunCollection: {self.run_name} with {len(self.samples)} read sets"

    def collect_samples(self, dir: Path) -> None:
        """Collects samples from the run directory"""
        fastq_files = list(dir.glob("*.fastq")) + list(dir.glob("*.fastq.gz"))
        for fq in fastq_files:
            sample_name, read_num, _ = extract_sample_info(fq)
            if sample_name not in self.samples:
                self.samples[sample_name] = RunSet(sample_name)
            if read_num == 1:
                self.samples[sample_name].read_1_fastqs.append(fq)
            elif read_num == 2:  #noqa
                self.samples[sample_name].read_2_fastqs.append(fq)

    def merge_all_samples(self) -> None:
        """Merges all lanes for each set of fastq files for all samples in the run collection"""
        for sample in self.samples.values():
            sample.merge_reads()


class RunSet:

    """Object to hold file and sample information for a single set of fastqs"""

    def __init__(self, sample_name: str) -> None:
        """Initializes a RunSet object"""
        self.sample_name = sample_name
        self.read_1_fastqs: list[Path] = []
        self.read_2_fastqs: list[Path] = []
        self.merged_r1: Path | None = None
        self.merged_r2: Path | None = None

    def __str__(self) -> str:
        """"""
        return self.sample_name

    def __repr__(self) -> str:
        """"""
        return f"RunSet: {self.sample_name} with {len(self.read_1_fastqs)} R1 fastq(s) and {len(self.read_2_fastqs)} R2 fastq(s)"

    def merge_reads(self) -> None:
        """Merges fastq files (lanes) for a given read number into a single file"""
        self.merged_r1 = self.read_1_fastqs[0].parent / f"{self.sample_name}_R1_merged.fastq.gz"
        self.merged_r2 = self.read_2_fastqs[0].parent / f"{self.sample_name}_R2_merged.fastq.gz"

        with Path(self.merged_r1).open("wb") as outfile:
            for fq in self.read_1_fastqs:
                with Path(fq).open("rb") as infile:
                    subprocess.run(["cat"], check=False, stdin=infile, stdout=outfile)

        with Path(self.merged_r2).open("wb") as outfile:
            for fq in self.read_2_fastqs:
                with Path(fq).open("rb") as infile:
                    subprocess.run(["cat"], check=False, stdin=infile, stdout=outfile)
