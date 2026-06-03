import csv
import gzip
import json
from pathlib import Path

from workflow_pathoscope.rust import PathoscopeResults, run_expectation_maximization


def _open_json(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")

    return open(path)


def _get_reference_otus(reference_data):
    if isinstance(reference_data, dict):
        return reference_data["otus"]

    return reference_data


def write_default_isolate_fasta(
    json_path: Path,
    target_path: Path,
) -> dict[str, int]:
    """Generate a FASTA file containing only default isolates from a reference JSON."""
    lengths = {}

    with _open_json(json_path) as f_json, open(target_path, "w") as f_target:
        for otu in _get_reference_otus(json.load(f_json)):
            for isolate in otu["isolates"]:
                if not isolate["default"]:
                    continue

                for sequence in isolate["sequences"]:
                    f_target.write(f">{sequence['_id']}\n{sequence['sequence']}\n")
                    lengths[sequence["_id"]] = len(sequence["sequence"])

    return lengths


def write_isolate_fasta(
    otu_ids: set[str],
    json_path: Path,
    target_path: Path,
) -> dict[str, int]:
    """Generate a FASTA file for all the isolates of the OTUs specified by ``otu_ids``.

    :param otu_ids: the list of OTU IDs for which to generate and index
    :param json_path: the path to the reference index json file
    :param target_path: the path to write the fasta file to
    :return: a dictionary of the lengths of all sequences keyed by their IDS

    """
    lengths = {}

    with _open_json(json_path) as f_json, open(target_path, "w") as f_target:
        for otu in _get_reference_otus(json.load(f_json)):
            if otu["_id"] in otu_ids:
                for isolate in otu["isolates"]:
                    for sequence in isolate["sequences"]:
                        f_target.write(f">{sequence['_id']}\n{sequence['sequence']}\n")
                        lengths[sequence["_id"]] = len(sequence["sequence"])

    return lengths


def write_report(
    path,
    pathoscope_results: PathoscopeResults,
):
    """Write pathoscope results to TSV report and return processed results dict.

    Float values are rounded to 10 decimal places for consistent precision
    across Rust/Python computations.
    """
    read_count = len(pathoscope_results.reads)

    tmp = zip(
        pathoscope_results.pi,
        pathoscope_results.refs,
        pathoscope_results.init_pi,
        pathoscope_results.best_hit_initial,
        pathoscope_results.best_hit_initial_reads,
        pathoscope_results.best_hit_final,
        pathoscope_results.best_hit_final_reads,
        pathoscope_results.level_1_initial,
        pathoscope_results.level_2_initial,
        pathoscope_results.level_1_final,
        pathoscope_results.level_2_final,
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
                len(pathoscope_results.refs),
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
                    "pi": round(x1[i], 10),
                    "best": round(x6[i], 10),
                    "high": round(x10[i], 10),
                    "low": round(x11[i], 10),
                    "reads": int(x7[i]),
                },
                "initial": {
                    "pi": round(x3[i], 10),
                    "best": round(x4[i], 10),
                    "high": round(x8[i], 10),
                    "low": round(x9[i], 10),
                    "reads": int(x5[i]),
                },
            }

    return results


def run_pathoscope(
    bam_path: Path,
    p_score_cutoff: float,
):
    """Run Pathoscope on an alignment file.

    Returns PathoscopeResults containing EM results and coverage data.

    :param alignment_path: The path to the SAM or BAM file.
    :param p_score_cutoff: The minimum allowed ``p_score`` for an alignment.
    """
    return run_expectation_maximization(
        str(bam_path),
        p_score_cutoff,
    )
