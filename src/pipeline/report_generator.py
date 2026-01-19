"""Report generation for pipeline results.

Generates developability scorecards and summary reports
for candidate CD3 binders.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import json
import datetime

from src.pipeline.filter_cascade import CandidateScore


@dataclass
class ReportConfig:
    """Configuration for report generation."""

    include_structures: bool = True
    include_sequences: bool = True
    include_filter_details: bool = True
    output_format: str = "html"  # "html", "json", or "both"


class ReportGenerator:
    """Generator for pipeline reports and scorecards."""

    def __init__(self, config: Optional[ReportConfig] = None):
        """Initialize report generator.

        Args:
            config: Report configuration.
        """
        self.config = config or ReportConfig()

    def generate_scorecard(
        self,
        candidate: CandidateScore,
        provenance: Optional[dict] = None,
    ) -> dict:
        """Generate developability scorecard for a candidate.

        Args:
            candidate: Candidate to score.
            provenance: Optional provenance metadata.

        Returns:
            Scorecard dictionary.
        """
        scorecard = candidate.to_dict()

        if provenance:
            scorecard["_provenance"] = provenance

        return scorecard

    def generate_summary_stats(
        self,
        candidates: list[CandidateScore],
        filter_stats: dict,
    ) -> dict:
        """Generate summary statistics for pipeline run.

        Args:
            candidates: Final candidates.
            filter_stats: Filtering statistics.

        Returns:
            Summary statistics dictionary.
        """
        # Epitope distribution
        epitope_counts = {}
        for c in candidates:
            epitope_counts[c.epitope_class] = epitope_counts.get(c.epitope_class, 0) + 1

        # Type distribution
        type_counts = {}
        for c in candidates:
            type_counts[c.binder_type] = type_counts.get(c.binder_type, 0) + 1

        # Source distribution
        source_counts = {}
        for c in candidates:
            source_counts[c.source] = source_counts.get(c.source, 0) + 1

        # Score statistics
        scores = [c.composite_score for c in candidates]
        pdockq_values = [c.pdockq for c in candidates if c.pdockq is not None]
        humanness_values = [c.oasis_score_mean or c.oasis_score_vh for c in candidates if c.oasis_score_mean or c.oasis_score_vh]

        return {
            "num_candidates": len(candidates),
            "filter_stats": filter_stats,
            "distributions": {
                "epitope": epitope_counts,
                "binder_type": type_counts,
                "source": source_counts,
            },
            "score_stats": {
                "composite": {
                    "min": min(scores) if scores else 0,
                    "max": max(scores) if scores else 0,
                    "mean": sum(scores) / len(scores) if scores else 0,
                },
                "pdockq": {
                    "min": min(pdockq_values) if pdockq_values else 0,
                    "max": max(pdockq_values) if pdockq_values else 0,
                    "mean": sum(pdockq_values) / len(pdockq_values) if pdockq_values else 0,
                },
                "humanness": {
                    "min": min(humanness_values) if humanness_values else 0,
                    "max": max(humanness_values) if humanness_values else 0,
                    "mean": sum(humanness_values) / len(humanness_values) if humanness_values else 0,
                },
            },
        }

    def generate_html_report(
        self,
        candidates: list[CandidateScore],
        filter_stats: dict,
        provenance: Optional[dict] = None,
        title: str = "CD3 Binder Design Report",
    ) -> str:
        """Generate HTML report.

        Args:
            candidates: Final candidates.
            filter_stats: Filtering statistics.
            provenance: Optional provenance metadata.
            title: Report title.

        Returns:
            HTML string.
        """
        summary = self.generate_summary_stats(candidates, filter_stats)

        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #333; }}
        h2 {{ color: #666; border-bottom: 1px solid #ccc; padding-bottom: 5px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .pass {{ color: green; }}
        .fail {{ color: red; }}
        .soft-fail {{ color: orange; }}
        .warning {{ background-color: #fff3cd; padding: 10px; border-radius: 5px; margin: 10px 0; }}
        .info {{ background-color: #cce5ff; padding: 10px; border-radius: 5px; margin: 10px 0; }}
        .metric {{ display: inline-block; margin: 10px; padding: 15px; background: #f5f5f5; border-radius: 5px; }}
        .metric-value {{ font-size: 24px; font-weight: bold; color: #333; }}
        .metric-label {{ font-size: 12px; color: #666; }}
    </style>
</head>
<body>
    <h1>{title}</h1>

    <div class="warning">
        <strong>Important:</strong> pDockQ is a structural confidence score, NOT an affinity predictor.
        All candidates require experimental validation.
    </div>

    <h2>Summary</h2>
    <div>
        <div class="metric">
            <div class="metric-value">{summary['num_candidates']}</div>
            <div class="metric-label">Final Candidates</div>
        </div>
        <div class="metric">
            <div class="metric-value">{filter_stats.get('total_input', 'N/A')}</div>
            <div class="metric-label">Total Designs</div>
        </div>
        <div class="metric">
            <div class="metric-value">{summary['score_stats']['pdockq']['mean']:.3f}</div>
            <div class="metric-label">Mean pDockQ</div>
        </div>
        <div class="metric">
            <div class="metric-value">{summary['score_stats']['humanness']['mean']:.3f}</div>
            <div class="metric-label">Mean Humanness</div>
        </div>
    </div>

    <h2>Candidate Rankings</h2>
    <table>
        <tr>
            <th>Rank</th>
            <th>ID</th>
            <th>Type</th>
            <th>Source</th>
            <th>pDockQ</th>
            <th>Humanness</th>
            <th>Epitope</th>
            <th>Score</th>
            <th>Flags</th>
        </tr>
"""

        for c in candidates:
            flags = ", ".join(c.risk_flags) if c.risk_flags else "None"
            pdockq_str = f"{c.pdockq:.3f}" if c.pdockq is not None else "N/A"
            humanness_str = f"{c.oasis_score_mean:.3f}" if c.oasis_score_mean is not None else "N/A"
            html += f"""        <tr>
            <td>{c.rank}</td>
            <td>{c.candidate_id}</td>
            <td>{c.binder_type}</td>
            <td>{c.source}</td>
            <td>{pdockq_str}</td>
            <td>{humanness_str}</td>
            <td>{c.epitope_class}</td>
            <td>{c.composite_score:.3f}</td>
            <td>{flags}</td>
        </tr>
"""

        html += """    </table>

    <h2>Distributions</h2>
    <h3>By Epitope Class</h3>
    <ul>
"""
        for epitope, count in summary['distributions']['epitope'].items():
            html += f"        <li>{epitope}: {count}</li>\n"

        html += """    </ul>

    <h3>By Binder Type</h3>
    <ul>
"""
        for btype, count in summary['distributions']['binder_type'].items():
            html += f"        <li>{btype}: {count}</li>\n"

        html += """    </ul>

    <h3>By Source</h3>
    <ul>
"""
        for source, count in summary['distributions']['source'].items():
            html += f"        <li>{source}: {count}</li>\n"

        html += """    </ul>

    <h2>Individual Scorecards</h2>
"""

        for c in candidates:
            scorecard = self.generate_scorecard(c, provenance)

            # Pre-format values that may be None
            pdockq_val = f"{c.pdockq:.3f}" if c.pdockq is not None else "N/A"
            humanness_val = f"{c.oasis_score_mean:.3f}" if c.oasis_score_mean is not None else "N/A"
            overlap_val = f"{c.okt3_overlap:.1%}" if c.okt3_overlap is not None else "N/A"
            binding_class = "pass" if c.filter_results.get("binding") else ""
            binding_status = c.filter_results.get("binding", "N/A") if hasattr(c, "filter_results") and c.filter_results else "N/A"
            humanness_status = c.filter_results.get("humanness", "N/A") if hasattr(c, "filter_results") and c.filter_results else "N/A"
            liabilities_status = c.filter_results.get("liabilities", "N/A") if hasattr(c, "filter_results") and c.filter_results else "N/A"
            seq_preview = c.sequence[:50] if c.sequence and len(c.sequence) > 50 else (c.sequence or "N/A")

            # Count liabilities safely
            deam_count = len(c.deamidation_sites) if c.deamidation_sites else 0
            glyc_count = len(c.glycosylation_sites) if c.glycosylation_sites else 0
            ox_count = len(c.oxidation_sites) if c.oxidation_sites else 0

            html += f"""    <h3>{c.candidate_id} (Rank {c.rank})</h3>
    <div class="info">
        <strong>Composite Score:</strong> {c.composite_score:.3f}<br>
        <strong>Sequence:</strong> <code>{seq_preview}...</code><br>
    </div>
    <table>
        <tr><th>Metric</th><th>Value</th><th>Status</th></tr>
        <tr>
            <td>pDockQ</td>
            <td>{pdockq_val}</td>
            <td class="{binding_class}">{binding_status}</td>
        </tr>
        <tr>
            <td>Humanness (OASis)</td>
            <td>{humanness_val}</td>
            <td>{humanness_status}</td>
        </tr>
        <tr>
            <td>Epitope Class</td>
            <td>{c.epitope_class}</td>
            <td>OKT3 overlap: {overlap_val}</td>
        </tr>
        <tr>
            <td>Liabilities</td>
            <td>Deam: {deam_count}, Glyc: {glyc_count}, Ox: {ox_count}</td>
            <td>{liabilities_status}</td>
        </tr>
    </table>
"""

        # Provenance
        if provenance:
            html += """    <h2>Provenance</h2>
    <ul>
"""
            for key, value in provenance.items():
                html += f"        <li><strong>{key}:</strong> {value}</li>\n"
            html += "    </ul>\n"

        html += f"""
    <hr>
    <p><em>Generated: {datetime.datetime.now().isoformat()}</em></p>
</body>
</html>
"""
        return html

    def generate_json_report(
        self,
        candidates: list[CandidateScore],
        filter_stats: dict,
        provenance: Optional[dict] = None,
    ) -> dict:
        """Generate JSON report.

        Args:
            candidates: Final candidates.
            filter_stats: Filtering statistics.
            provenance: Optional provenance metadata.

        Returns:
            Report dictionary.
        """
        summary = self.generate_summary_stats(candidates, filter_stats)

        return {
            "summary": summary,
            "candidates": [self.generate_scorecard(c, provenance) for c in candidates],
            "provenance": provenance,
            "generated_at": datetime.datetime.now().isoformat(),
        }

    def save_report(
        self,
        candidates: list[CandidateScore],
        filter_stats: dict,
        output_dir: str,
        provenance: Optional[dict] = None,
        timestamp: Optional[str] = None,
    ) -> list[str]:
        """Save reports to files.

        Args:
            candidates: Final candidates.
            filter_stats: Filtering statistics.
            output_dir: Output directory.
            provenance: Optional provenance metadata.
            timestamp: Optional timestamp for filename.

        Returns:
            List of saved file paths.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        timestamp = timestamp or datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        saved_files = []

        # JSON report
        if self.config.output_format in ["json", "both"]:
            json_report = self.generate_json_report(candidates, filter_stats, provenance)
            json_path = output_dir / f"report_{timestamp}.json"
            with open(json_path, "w") as f:
                json.dump(json_report, f, indent=2)
            saved_files.append(str(json_path))

        # HTML report
        if self.config.output_format in ["html", "both"]:
            html_report = self.generate_html_report(candidates, filter_stats, provenance)
            html_path = output_dir / f"report_{timestamp}.html"
            with open(html_path, "w") as f:
                f.write(html_report)
            saved_files.append(str(html_path))

        # Individual scorecards
        scorecards_dir = output_dir / "scorecards"
        scorecards_dir.mkdir(exist_ok=True)

        for candidate in candidates:
            scorecard = self.generate_scorecard(candidate, provenance)
            scorecard_path = scorecards_dir / f"{candidate.candidate_id}_scorecard.json"
            with open(scorecard_path, "w") as f:
                json.dump(scorecard, f, indent=2)
            saved_files.append(str(scorecard_path))

        return saved_files


def generate_report(
    candidates: list[CandidateScore],
    filter_stats: dict,
    output_dir: str,
    provenance: Optional[dict] = None,
) -> list[str]:
    """Convenience function for report generation.

    Args:
        candidates: Final candidates.
        filter_stats: Filtering statistics.
        output_dir: Output directory.
        provenance: Optional provenance metadata.

    Returns:
        List of saved file paths.
    """
    config = ReportConfig(output_format="both")
    generator = ReportGenerator(config)

    saved = generator.save_report(
        candidates=candidates,
        filter_stats=filter_stats,
        output_dir=output_dir,
        provenance=provenance,
    )

    print(f"Reports saved: {len(saved)} files")
    for path in saved[:3]:  # Show first 3
        print(f"  {path}")
    if len(saved) > 3:
        print(f"  ... and {len(saved) - 3} more")

    return saved
