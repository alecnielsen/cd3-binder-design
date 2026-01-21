from dataclasses import dataclass

from src.analysis.liabilities import LiabilityReport, LiabilitySite
from src.analysis.humanness import HumannessReport, PairedHumannessReport
from src.pipeline.design_pipeline import DesignPipeline


@dataclass
class _DummyPhysicochemical:
    net_charge: float = 0.0
    isoelectric_point: float = 7.0


@dataclass
class _DummyAggregation:
    hydrophobic_patches: int = 0


def test_scfv_liability_offsets_include_linker(monkeypatch):
    pipeline = DesignPipeline()

    def fake_scan_with_cdr_detection(self, sequence, chain_type="H", numbering_scheme="imgt"):
        if chain_type == "H":
            self.cdr_positions = {"H3": (1, 2)}
            return LiabilityReport(
                sequence=sequence,
                deamidation_sites=[
                    LiabilitySite(
                        motif="NG",
                        position=0,
                        liability_type="deamidation",
                        in_cdr=True,
                        cdr_name="H3",
                    )
                ],
                isomerization_sites=[],
                glycosylation_sites=[],
                oxidation_sites=[],
                unpaired_cysteines=0,
                total_liabilities=1,
                cdr_liabilities=1,
                cdr_oxidation_count=0,
            )
        self.cdr_positions = {"L3": (1, 2)}
        return LiabilityReport(
            sequence=sequence,
            deamidation_sites=[
                LiabilitySite(
                    motif="NS",
                    position=1,
                    liability_type="deamidation",
                    in_cdr=True,
                    cdr_name="L3",
                )
            ],
            isomerization_sites=[],
            glycosylation_sites=[],
            oxidation_sites=[],
            unpaired_cysteines=0,
            total_liabilities=1,
            cdr_liabilities=1,
            cdr_oxidation_count=0,
        )

    def fake_assess(self, vh_sequence, vl_sequence=None, **_kwargs):
        return type(
            "DummyDevReport",
            (),
            {
                "cdr_h3_length": 2,
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

    vh = "AAAA"
    vl = "BBBB"
    candidates = [{"vh": vh, "vl": vl, "binder_type": "scfv"}]
    results = pipeline.run_analysis(candidates)
    score = results[0]

    linker_len = len(pipeline.config.formatting.scfv_linker)
    vl_offset = len(vh) + linker_len

    assert score.full_sequence == vh + pipeline.config.formatting.scfv_linker + vl
    assert score.deamidation_sites == [0, 1 + vl_offset]
    assert score.cdr_positions["L3"] == (1 + vl_offset, 2 + vl_offset)
