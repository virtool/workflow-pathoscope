import csv
from functools import cached_property
from pathlib import Path
from typing import Any, Generator

import rust_utils


class SamLine:
    def __init__(self, line: str):
        self._line = line

    def __str__(self) -> str:
        return self.line

    @property
    def line(self) -> str:
        """The SAM line used to create the object."""
        return self._line

    @property
    def read_id(self) -> str:
        """The ID of the mapped read."""
        return self.fields[0]

    @cached_property
    def read_length(self) -> int:
        """The length of the mapped read."""
        return len(self.fields[9])

    @cached_property
    def fields(self) -> list[Any]:
        """The SAM fields"""
        return self.line.split("\t")

    @cached_property
    def position(self) -> int:
        """The position of the read on the reference."""
        return int(self.fields[3])

    @cached_property
    def score(self) -> float:
        """The Pathoscope score for the alignment."""
        return find_sam_align_score(self.fields)

    @cached_property
    def bitwise_flag(self) -> int:
        """The SAM bitwise flag."""
        return int(self.fields[1])

    @cached_property
    def unmapped(self) -> bool:
        """The read is unmapped.

        This value is derived from the bitwise flag (0x4: segment unmapped).
        """
        return self.bitwise_flag & 4 == 4

    @cached_property
    def ref_id(self) -> str:
        """The ID of the mapped reference sequence."""
        return self.fields[2]


def calculate_coverage(sam_path: Path, ref_lengths: list[str, int]):
    coverage_dict = {}
    pos_length_list = []

    for line in parse_sam(sam_path):
        coverage_dict[line.ref_id] = [0] * ref_lengths[line.ref_id]
        pos_length_list.append((line.ref_id, line.position, line.read_length))

    for ref_id, pos, length in pos_length_list:
        start_index = pos - 1

        for i in range(start_index, start_index + length):
            try:
                coverage_dict[ref_id][i] += 1
            except IndexError:
                pass

    return coverage_dict


def find_sam_align_score(fields: list[Any]) -> float:
    """Find the Bowtie2 alignment score for the given split line (``fields``).

    Searches the SAM fields for the ``AS:i`` substring and extracts the Bowtie2-specific
    alignment score. This will not work for other aligners.

    :param fields: a SAM line that has been split on "\t"
    :return: the alignment score

    """
    read_length = float(len(fields[9]))

    for field in fields:
        if field.startswith("AS:i:"):
            a_score = int(field[5:])
            return a_score + read_length

    raise ValueError("Could not find alignment score")


def parse_sam(
    path: Path,
    p_score_cutoff: float = 0.01,
) -> Generator[SamLine, None, None]:
    """Parse a SAM file and yield :class:`SamLine` objects.

    :param path: The path to the SAM file.
    :param p_score_cutoff: The minimum allowed ``p_score`` for an alignment.
    :return: A generator of sam lines.

    """
    with open(path) as f:
        for line in f:
            if line[0] == "#" or line[0] == "@":
                continue

            sam_line = SamLine(line)

            if sam_line.unmapped:
                continue

            if sam_line.score < p_score_cutoff:
                continue

            yield SamLine(line)


def write_report(
    path,
    pi,
    refs,
    read_count,
    init_pi,
    best_hit_initial,
    best_hit_initial_reads,
    best_hit_final,
    best_hit_final_reads,
    level_1_initial,
    level_2_initial,
    level_1_final,
    level_2_final,
):
    tmp = zip(
        pi,
        refs,
        init_pi,
        best_hit_initial,
        best_hit_initial_reads,
        best_hit_final,
        best_hit_final_reads,
        level_1_initial,
        level_2_initial,
        level_1_final,
        level_2_final,
    )

    tmp = sorted(tmp, reverse=True)

    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 = zip(*tmp)

    end = 0

    for i, _ in enumerate(x10):
        if x1[i] < 0.01 and x10[i] <= 0 and x11[i] <= 0:
            break

        end += 1

    with open(path, "w") as handle:
        csv_writer = csv.writer(handle, delimiter="\t")

        csv_writer.writerow(
            [
                "Total Number of Aligned Reads:",
                read_count,
                "Total Number of Mapped Genomes:",
                len(refs),
            ],
        )

        csv_writer.writerow(
            [
                "Genome",
                "Final Guess",
                "Final Best Hit",
                "Final Best Hit Read Numbers",
                "Final High Confidence Hits",
                "Final Low Confidence Hits",
                "Initial Guess",
                "Initial Best Hit",
                "Initial Best Hit Read Numbers",
                "Initial High Confidence Hits",
                "Initial Low Confidence Hits",
            ],
        )

        # Change the column order with zip.
        csv_writer.writerows(
            zip(
                x2[:end],
                x1[:end],
                x6[:end],
                x7[:end],
                x10[:end],
                x11[:end],
                x3[:end],
                x4[:end],
                x5[:end],
                x8[:end],
                x9[:end],
            ),
        )

    results = {}

    for i, ref_id in enumerate(x2[:end]):
        if x1[i] < 0.01 and x10[i] <= 0 and x11[i] <= 0:
            pass
        else:
            results[ref_id] = {
                "final": {
                    "pi": x1[i],
                    "best": x6[i],
                    "high": x10[i],
                    "low": x11[i],
                    "reads": int(x7[i]),
                },
                "initial": {
                    "pi": x3[i],
                    "best": x4[i],
                    "high": x8[i],
                    "low": x9[i],
                    "reads": int(x5[i]),
                },
            }

    return results


def run_pathoscope(sam_path: Path, reassigned_path: Path, p_score_cutoff: float):
    """Run Pathoscope on the SAM file at ``sam_path`` with the given ``p_score_cutoff``.

    The output is a reassigned SAM file at ``reassigned_path``. The SAM file can be
    used to generate a report with :func:`write_report`.

    :param sam_path: The path to the SAM file.
    :param reassigned_path: The path to the reassigned SAM file.
    :param p_score_cutoff: The minimum allowed ``p_score`` for an alignment.
    """
    return rust_utils.run_expectation_maximization(
        str(sam_path),
        str(reassigned_path),
        p_score_cutoff,
    )
