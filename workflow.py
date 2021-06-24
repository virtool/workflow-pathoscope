from pathlib import Path
from types import SimpleNamespace
from typing import List

import aiofiles
from virtool_workflow import fixture, step
from virtool_workflow.analysis.indexes import Index
from virtool_workflow.analysis.reads import Reads
from virtool_workflow.execution.run_subprocess import RunSubprocess


@fixture
def index(indexes: List[Index]):
    return indexes[0]


@fixture
def intermediate():
    """A namespace for storing intermediate values."""
    return SimpleNamespace(
        to_otus=set(),
        lengths=None,
    )


@fixture
def isolate_path(work_path: Path):
    path = work_path / "isolates"
    path.mkdir()
    return path


@fixture
def p_score_cutoff():
    return 0.01


@step
async def map_default_isolates(
        pathoscope,
        intermediate: dict,
        reads: Reads,
        index: Index,
        proc: int,
        p_score_cutoff: float,
        run_subprocess: RunSubprocess,
):
    """
    Map reads to the main OTU reference.

    This will be used to identify canididate OTUs.
    """

    async def stdout_handler(line: str):
        line = line.decode()

        if line[0] == "#" or line[0] == "@":
            return

        fields = line.split("\t")

        # Bitwise FLAG - 0x4: segment unmapped
        if int(fields[1]) & 0x4 == 4:
            return

        ref_id = fields[2]

        if ref_id == "*":
            return

        # Skip if the p_score does not meet the minimum cutoff.
        if pathoscope.find_sam_align_score(fields) < p_score_cutoff:
            return

        intermediate.to_otus.add(ref_id)

    await run_subprocess(
        [
            "bowtie2",
            "-p", str(proc),
            "--no-unal"
            "--local",
            "--score-min", "L,20,1.0",
            "-N", "0",
            "-L", "15",
            "-x", index.path,
            "-U", f"{reads.left},{reads.right}",
        ],
        wait=True,
        stdout_handler=stdout_handler
    )

    return f"Mapped reats to OTUs {intermediate.to_otus}"


@step
async def write_isolate_fasta(
    intermediate,
    index: Index,
    proc: int,
    isolate_path: Path,
):
    fasta_path = isolate_path/"isolate_index.fa",
    intermediate.lengths = await index.write_isolate_fasta(
        [index.get_otu_id_by_sequence_id(id_) for id_ in intermediate.to_otus],
        fasta_path,
        proc
    )

    intermediate.isolate_fasta_path = fasta_path

    return "Produced isolate fasta file."

