import copy
import csv
import math
import rust_utils
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, Generator, List


class SamLine:
    def __init__(self, line: str):
        self._line = line

    def __str__(self) -> str:
        return self.line

    @property
    def line(self) -> str:
        """
        The SAM line used to create the object.
        """
        return self._line

    @property
    def read_id(self) -> str:
        """
        The ID of the mapped read.
        """
        return self.fields[0]

    @cached_property
    def read_length(self) -> int:
        """
        The length of the mapped read.
        """
        return len(self.fields[9])

    @cached_property
    def fields(self) -> List[Any]:
        """
        The SAM fields
        """
        return self.line.split("\t")

    @cached_property
    def position(self) -> int:
        """
        The position of the read on the reference.
        """
        return int(self.fields[3])

    @cached_property
    def score(self) -> float:
        """
        The Pathoscope score for the alignment.
        """
        return find_sam_align_score(self.fields)

    @cached_property
    def bitwise_flag(self) -> int:
        """
        The SAM bitwise flag.
        """
        return int(self.fields[1])

    @cached_property
    def unmapped(self) -> bool:
        """
        The read is unmapped.

        This value is derived from the bitwise flag (0x4: segment unmapped).
        """
        return self.bitwise_flag & 4 == 4

    @cached_property
    def ref_id(self) -> str:
        """
        The ID of the mapped reference sequence.
        """
        return self.fields[2]


def parse_sam(
    path: Path, p_score_cutoff: float = 0.01
) -> Generator[SamLine, None, None]:
    """
    Parse a SAM file and yield :class:`SamLine` objects.

    :param path: The path to the SAM file.
    :param p_score_cutoff: The minimum allowed ``p_score`` for an alignment.
    :return: A generator of sam lines.

    """
    with open(path, "r") as f:
        for line in f:
            if line[0] == "#" or line[0] == "@":
                continue

            sam_line = SamLine(line)

            if sam_line.unmapped:
                continue

            if sam_line.score < p_score_cutoff:
                continue

            yield SamLine(line)


def rescale_samscore(u, nu, max_score, min_score):
    if min_score < 0:
        scaling_factor = 100.0 / max_score - min_score
    else:
        scaling_factor = 100.0 / max_score

    for read_index in u:
        if min_score < 0:
            u[read_index][1][0] = u[read_index][1][0] - min_score

        u[read_index][1][0] = math.exp(u[read_index][1][0] * scaling_factor)
        u[read_index][3] = u[read_index][1][0]

    for read_index in nu:
        nu[read_index][3] = 0.0

        for i in range(0, len(nu[read_index][1])):
            if min_score < 0:
                nu[read_index][1][i] = nu[read_index][1][i] - min_score

            nu[read_index][1][i] = math.exp(nu[read_index][1][i] * scaling_factor)

            if nu[read_index][1][i] > nu[read_index][3]:
                nu[read_index][3] = nu[read_index][1][i]

    return u, nu


def find_sam_align_score(fields: List[Any]) -> float:
    """
    Find the Bowtie2 alignment score for the given split line (``fields``).

    Searches the SAM fields for the ``AS:i`` substring and extracts the Bowtie2-specific alignment score. This will not
    work for other aligners.

    :param fields: a SAM line that has been split on "\t"
    :return: the alignment score

    """
    read_length = float(len(fields[9]))

    for field in fields:
        if field.startswith("AS:i:"):
            a_score = int(field[5:])
            return a_score + read_length

    raise ValueError("Could not find alignment score")


def build_matrix(sam_path: Path, p_score_cutoff=0.01):
    u = dict()
    nu = dict()

    h_read_id = {}
    h_ref_id = {}

    refs = []
    reads = []

    ref_count = 0
    read_count = 0

    max_score = 0
    min_score = 0

    for sam_line in parse_sam(sam_path, p_score_cutoff):
        if sam_line.score < p_score_cutoff:
            continue

        min_score = min(min_score, sam_line.score)
        max_score = max(max_score, sam_line.score)

        ref_index = h_ref_id.get(sam_line.ref_id, -1)

        if ref_index == -1:
            ref_index = ref_count
            h_ref_id[sam_line.ref_id] = ref_index
            refs.append(sam_line.ref_id)
            ref_count += 1

        read_index = h_read_id.get(sam_line.read_id, -1)

        if read_index == -1:
            # hold on this new read. first, wrap previous read profile and see if any
            # previous read has a same profile with that!
            read_index = read_count
            h_read_id[sam_line.read_id] = read_index
            reads.append(sam_line.read_id)
            read_count += 1
            u[read_index] = [
                [ref_index],
                [sam_line.score],
                [float(sam_line.score)],
                sam_line.score,
            ]
        else:
            if read_index in u:
                if ref_index in u[read_index][0]:
                    continue
                nu[read_index] = u[read_index]
                del u[read_index]

            if ref_index in nu[read_index][0]:
                continue

            nu[read_index][0].append(ref_index)
            nu[read_index][1].append(sam_line.score)

            if sam_line.score > nu[read_index][3]:
                nu[read_index][3] = sam_line.score

    u, nu = rescale_samscore(u, nu, max_score, min_score)

    for read_index in u:
        # keep ref_index and score only
        u[read_index] = [u[read_index][0][0], u[read_index][1][0]]

    for read_index in nu:
        p_score_sum = sum(nu[read_index][1])
        # Normalize p_score.
        nu[read_index][2] = [k / p_score_sum for k in nu[read_index][1]]

    return u, nu, refs, reads


def em(u, nu, genomes, max_iter, epsilon, pi_prior, theta_prior):
    genome_count = len(genomes)

    pi = [1.0 / genome_count] * genome_count
    init_pi = copy.copy(pi)
    theta = copy.copy(pi)

    pi_sum_0 = [0] * genome_count

    u_weights = [u[i][1] for i in u]

    max_u_weights = 0
    u_total = 0

    if u_weights:
        max_u_weights = max(u_weights)
        u_total = sum(u_weights)

    for i in u:
        pi_sum_0[u[i][0]] += u[i][1]

    nu_weights = [nu[i][3] for i in nu]

    max_nu_weights = 0
    nu_total = 0

    if nu_weights:
        max_nu_weights = max(nu_weights)
        nu_total = sum(nu_weights)

    prior_weight = max(max_u_weights, max_nu_weights)
    nu_length = len(nu)

    if nu_length == 0:
        nu_length = 1

    # EM iterations
    for i in range(max_iter):
        pi_old = pi
        theta_sum = [0 for _ in genomes]

        # E Step
        for j in nu:
            z = nu[j]

            # A set of any genome mapping with j
            ind = z[0]

            # Get relevant pis for the read
            pi_tmp = [pi[k] for k in ind]

            # Get relevant thetas for the read.
            theta_tmp = [theta[k] for k in ind]

            # Calculate non-normalized xs
            x_tmp = [1.0 * pi_tmp[k] * theta_tmp[k] * z[1][k] for k in range(len(ind))]

            x_sum = sum(x_tmp)

            # Avoid dividing by 0 at all times.
            if x_sum == 0:
                x_norm = [0.0 for _ in x_tmp]
            else:
                # Normalize new xs.
                x_norm = [1.0 * k / x_sum for k in x_tmp]

            # Update x in nu.
            nu[j][2] = x_norm

            for k, _ in enumerate(ind):
                # Keep weighted running tally for theta
                theta_sum[ind[k]] += x_norm[k] * nu[j][3]

        # M step
        pi_sum = [theta_sum[k] + pi_sum_0[k] for k in range(len(theta_sum))]
        pip = pi_prior * prior_weight

        # Update pi.
        pi = [
            (1.0 * k + pip) / (u_total + nu_total + pip * len(pi_sum)) for k in pi_sum
        ]

        if i == 0:
            init_pi = pi

        theta_p = theta_prior * prior_weight

        nu_total_div = nu_total

        if nu_total_div == 0:
            nu_total_div = 1

        theta = [
            (1.0 * k + theta_p) / (nu_total_div + theta_p * len(theta_sum))
            for k in theta_sum
        ]

        cutoff = 0.0

        for k, _ in enumerate(pi):
            cutoff += abs(pi_old[k] - pi[k])

        if cutoff <= epsilon or nu_length == 1:
            break

    return init_pi, pi, theta, nu


def compute_best_hit(u, nu, refs, reads):
    ref_count = len(refs)

    best_hit_reads = [0.0] * ref_count
    level_1_reads = [0.0] * ref_count
    level_2_reads = [0.0] * ref_count

    for i in u:
        best_hit_reads[u[i][0]] += 1
        level_1_reads[u[i][0]] += 1

    for j in nu:
        z = nu[j]
        ind = z[0]
        x_norm = z[2]
        best_ref = max(x_norm)
        num_best_ref = 0

        for i, _ in enumerate(x_norm):
            if x_norm[i] == best_ref:
                num_best_ref += 1

        num_best_ref = num_best_ref or 1

        for i, _ in enumerate(x_norm):
            if x_norm[i] == best_ref:
                best_hit_reads[ind[i]] += 1.0 / num_best_ref

                if x_norm[i] >= 0.5:
                    level_1_reads[ind[i]] += 1
                elif x_norm[i] >= 0.01:
                    level_2_reads[ind[i]] += 1

    ref_count = len(refs)
    read_count = len(reads)

    best_hit = [best_hit_reads[k] / read_count for k in range(ref_count)]
    level_1 = [level_1_reads[k] / read_count for k in range(ref_count)]
    level_2 = [level_2_reads[k] / read_count for k in range(ref_count)]

    return best_hit_reads, best_hit, level_1, level_2


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

    no_cutoff = False

    for i, _ in enumerate(x10):
        if not no_cutoff and x1[i] < 0.01 and x10[i] <= 0 and x11[i] <= 0:
            break

        if i == (len(x10) - 1):
            i += 1

    # Changing the column order here
    tmp = zip(
        x2[:i],
        x1[:i],
        x6[:i],
        x7[:i],
        x10[:i],
        x11[:i],
        x3[:i],
        x4[:i],
        x5[:i],
        x8[:i],
        x9[:i],
    )

    with open(path, "w") as handle:
        csv_writer = csv.writer(handle, delimiter="\t")

        header = [
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
        ]

        header1 = [
            "Total Number of Aligned Reads:",
            read_count,
            "Total Number of Mapped Genomes:",
            len(refs),
        ]

        csv_writer.writerow(header1)
        csv_writer.writerow(header)
        csv_writer.writerows(tmp)

    results = dict()

    for i, ref_id in enumerate(x2[:i]):
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


def calculate_coverage(sam_path: Path, ref_lengths: Dict[str, int]):
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


def run(sam_path: Path, reassigned_path: Path, p_score_cutoff: float):
    # rust binding to improve runtime performance
    # computes:
    #   buildMatrix
    #   em
    #   rewriteAlign
    #   computeBestHit
    #   + adjacent code
    return rust_utils.run_expectation_maximization(
        str(sam_path), str(reassigned_path), p_score_cutoff
    )
