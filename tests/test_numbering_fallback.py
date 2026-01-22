from src.analysis import numbering


def test_get_cdr_positions_returns_empty_on_numbering_failure(monkeypatch):
    monkeypatch.setattr(numbering, "number_sequence", lambda *args, **_kwargs: None)

    positions = numbering.get_cdr_positions("AAAA", chain_type="H")

    assert positions == {}
