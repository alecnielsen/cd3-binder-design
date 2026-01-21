from src.pipeline.filter_cascade import CandidateScore, FilterCascade, FilterResult


def test_filter_aggregation_uses_cdr_positions():
    cascade = FilterCascade()
    candidate = CandidateScore(
        candidate_id="cdr_aromatic",
        sequence="FWYFWAAAAAAAAAAAAAAAAAAAA",
    )
    candidate.full_sequence = candidate.sequence
    candidate.cdr_positions = {"H1": (0, 4)}

    result = cascade.filter_aggregation(candidate)

    assert result == FilterResult.SOFT_FAIL


def test_filter_aggregation_ignores_framework_aromatics_when_cdrs_clean():
    cascade = FilterCascade()
    candidate = CandidateScore(
        candidate_id="framework_aromatic",
        sequence="FWYFWFWYFW" + "AAAAA",
    )
    candidate.full_sequence = candidate.sequence
    # Define a clean CDR region with no aromatics.
    candidate.cdr_positions = {"H1": (10, 14)}

    result = cascade.filter_aggregation(candidate)

    assert result == FilterResult.PASS
