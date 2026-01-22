from src.pipeline.filter_cascade import CandidateScore, FilterCascade, FilterResult


def test_filter_humanness_treats_zero_as_real_score():
    cascade = FilterCascade()
    candidate = CandidateScore(
        candidate_id="zero_humanness",
        sequence="AAAA",
    )
    candidate.oasis_score_mean = 0.0

    result = cascade.filter_humanness(candidate)

    assert result == FilterResult.FAIL
