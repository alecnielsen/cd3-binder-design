from dataclasses import dataclass

from src.analysis.liabilities import LiabilityReport
from src.analysis.humanness import HumannessReport, PairedHumannessReport
from src.pipeline.design_pipeline import DesignPipeline


@dataclass
class _DummyPhysicochemical:
    net_charge: float = 0.0
    isoelectric_point: float = 7.0


@dataclass
class _DummyAggregation:
    hydrophobic_patches: int = 0


def test_run_analysis_parses_scfv_sequence_when_binder_type_set(monkeypatch):
    pipeline = DesignPipeline()

    def fake_scan_with_cdr_detection(self, sequence, chain_type="H", numbering_scheme="imgt"):
        self.cdr_positions = {}
        return LiabilityReport(
            sequence=sequence,
            deamidation_sites=[],
            isomerization_sites=[],
            glycosylation_sites=[],
            oxidation_sites=[],
            unpaired_cysteines=0,
            total_liabilities=0,
            cdr_liabilities=0,
            cdr_oxidation_count=0,
        )

    def fake_assess(self, vh_sequence, vl_sequence=None, **_kwargs):
        return type(
            "DummyDevReport",
            (),
            {
                "cdr_h3_length": 10,
                "physicochemical": _DummyPhysicochemical(),
                "aggregation": _DummyAggregation(),
            },
        )()

    def fake_humanness_pair(vh_sequence, vl_sequence=None):
        return PairedHumannessReport(
            vh_report=HumannessReport(
                sequence=vh_sequence,
                chain_type="H",
                oasis_score=0.9,
            ),
            vl_report=HumannessReport(
                sequence=vl_sequence or "",
                chain_type="L",
                oasis_score=0.9,
            )
            if vl_sequence
            else None,
            mean_score=0.9,
        )

    monkeypatch.setattr(
        "src.analysis.liabilities.LiabilityScanner.scan_with_cdr_detection",
        fake_scan_with_cdr_detection,
    )
    monkeypatch.setattr(
        "src.analysis.developability.DevelopabilityAssessor.assess",
        fake_assess,
    )
    monkeypatch.setattr(
        "src.analysis.humanness.score_humanness_pair",
        fake_humanness_pair,
    )

    vh = "A" * 120
    vl = "B" * 110
    linker = pipeline.config.formatting.scfv_linker
    scfv = vh + linker + vl

    candidates = [{"sequence": scfv, "binder_type": "scfv"}]
    results = pipeline.run_analysis(candidates)
    score = results[0]

    assert score.sequence == vh
    assert score.sequence_vl == vl
    assert score.full_sequence == scfv
