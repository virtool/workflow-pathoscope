import pickle
from pathlib import Path

import pytest

import pathoscope

BASE_PATH = Path.cwd() / "test_files"
BEST_HIT_PATH = Path.cwd() / BASE_PATH / "best_hit"
EM_PATH = BASE_PATH / "em"
MATRIX_PATH = BASE_PATH / "ps_matrix"
SAM_PATH = BASE_PATH / "test_al.sam"
SCORES = BASE_PATH / "scores"
TSV_PATH = BASE_PATH / "report.tsv"
UNU_PATH = BASE_PATH / "unu"
VTA_PATH = BASE_PATH / "test.vta"


@pytest.fixture(scope="session")
def expected_em():
    with open(EM_PATH, "rb") as handle:
        em_dict = pickle.load(handle)

    return em_dict


@pytest.fixture(scope="session")
def expected_scores():
    with open(SCORES, "rb") as handle:
        scores = pickle.load(handle)

    return scores


def test_find_sam_align_score(sam_line, expected_scores):
    assert (
        pathoscope.find_sam_align_score(sam_line) == expected_scores["".join(sam_line)]
    )


@pytest.mark.datafiles(SAM_PATH)
def test_build_matrix(data_regression, datafiles):
    sam_path = datafiles / "test_al.sam"

    with open(MATRIX_PATH, "rb") as handle:
        expected = pickle.load(handle)

    actual = pathoscope.build_matrix(sam_path, 0.01)
    u, nu, refs, reads = actual

    data_regression.check([u, nu, sorted(refs), sorted(reads)])

    assert expected[0] == actual[0]
    assert expected[1] == actual[1]
    assert sorted(expected[2]) == sorted(actual[2])
    assert sorted(expected[3]) == sorted(actual[3])


@pytest.mark.datafiles(SAM_PATH)
@pytest.mark.parametrize("theta_prior", [0])
@pytest.mark.parametrize("pi_prior", [0])
@pytest.mark.parametrize("epsilon", [1e-6])
@pytest.mark.parametrize("max_iter", [5])
def test_em(
    data_regression,
    datafiles,
    tmp_path,
    theta_prior,
    pi_prior,
    epsilon,
    max_iter,
    expected_em,
):
    sam_path = datafiles / "test_al.sam"

    u, nu, refs, _ = pathoscope.build_matrix(sam_path, 0.01)

    result = pathoscope.em(u, nu, refs, max_iter, epsilon, pi_prior, theta_prior)
    init_pi, pi, theta, updated_nu = result

    file_string = "_".join(
        [str(i) for i in ["em", theta_prior, pi_prior, epsilon, max_iter]]
    )

    for i in [0, 1, 2]:
        assert sorted(result[i]) == sorted(expected_em[file_string][i])

    assert result[3] == expected_em[file_string][3]

    data_regression.check([sorted(init_pi), sorted(pi), sorted(theta), updated_nu])


def test_compute_best_hit(data_regression):
    """
    Test that :meth:`compute_best_hit` gives the expected result given some input data.

    """
    with open(MATRIX_PATH, "rb") as handle:
        matrix_tuple = pickle.load(handle)

    data_regression.check(pathoscope.compute_best_hit(*matrix_tuple))


@pytest.mark.datafiles(SAM_PATH)
def test_rewrite_align(datafiles, file_regression, tmp_path):
    with open(UNU_PATH, "rb") as f:
        u, nu = pickle.load(f)

    rewrite_path = tmp_path / "rewrite.sam"

    pathoscope.rewrite_align(u, nu, datafiles / "test_al.sam", 0.01, rewrite_path)

    with open(rewrite_path, "r") as f:
        file_regression.check(f.read())


@pytest.mark.datafiles(SAM_PATH)
def test_calculate_coverage(datafiles, tmp_path, sam_path):
    ref_lengths = dict()

    sam_path = datafiles / "test_al.sam"

    with open(sam_path, "r") as handle:
        for line in handle:
            if line[0:3] == "@SQ":
                ref_id = None
                length = None

                for field in line.rstrip().split("\t"):
                    if "SN:" in field:
                        ref_id = field.split(":")[1]
                    if "LN:" in field:
                        length = int(field.split(":")[1])

                assert ref_id and length
                assert ref_id not in ref_lengths

                ref_lengths[ref_id] = length

    pathoscope.calculate_coverage(sam_path, ref_lengths)


@pytest.mark.datafiles(SAM_PATH)
def test_write_report(datafiles, file_regression, tmp_path):
    sam_path = datafiles / "test_al.sam"

    u, nu, refs, reads = pathoscope.build_matrix(sam_path, 0.01)

    with open(BEST_HIT_PATH, "rb") as handle:
        (
            best_hit_initial_reads,
            best_hit_initial,
            level_1_initial,
            level_2_initial,
        ) = pickle.load(handle)

    init_pi, pi, _, nu = pathoscope.em(u, nu, refs, 30, 1e-7, 0, 0)

    (
        best_hit_final_reads,
        best_hit_final,
        level_1_final,
        level_2_final,
    ) = pathoscope.compute_best_hit(u, nu, refs, reads)

    report_path = tmp_path / "report.tsv"

    pathoscope.write_report(
        report_path,
        pi,
        refs,
        len(reads),
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

    with open(report_path, "r") as f:
        file_regression.check(f.read())
