"""Asymmetric Fab + VHH bispecific format assembly.

One arm is a standard Fab (for tumor target), the other arm is a VHH nanobody
(for CD3) fused directly to the Fc. Uses knob-in-hole for heterodimerization.

Structure:
- Heavy chain 1 (Knob): VH(target) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2 (Hole): VHH(CD3) - Hinge - CH2 - CH3(Hole)
- Light chain: VL(target) - CL

Advantages:
- Simplest asymmetric format (single-domain CD3 binder)
- Smallest molecular weight
- No light chain pairing issues for CD3 arm
"""

from typing import Optional

from src.formatting.base import (
    BispecificFormatter,
    BispecificConstruct,
    AntibodyChain,
    SequenceLibrary,
)


class FabVhhFormatter(BispecificFormatter):
    """Formatter for asymmetric Fab + VHH bispecifics."""

    @property
    def format_name(self) -> str:
        return "fab_vhh"

    def assemble(
        self,
        target_vh: str,
        target_vl: str,
        cd3_binder: str,  # VHH for CD3
        cd3_binder_vl: Optional[str] = None,  # Ignored for VHH
        name: str = "fab_vhh_bispecific",
        target_name: str = "HER2",
    ) -> BispecificConstruct:
        """Assemble Fab + VHH bispecific.

        Args:
            target_vh: VH sequence for tumor target Fab arm.
            target_vl: VL sequence for tumor target Fab arm.
            cd3_binder: VHH sequence for CD3 arm.
            cd3_binder_vl: Ignored (VHH has no VL).
            name: Name for the construct.
            target_name: Name of tumor target.

        Returns:
            BispecificConstruct with 3 chains.
        """
        # Get constant region sequences
        ch1 = self.library.get_ch1()
        hinge = self.library.get_hinge()
        ch2 = self.library.get_ch2()
        ch3_knob = self.library.get_ch3_knob()
        ch3_hole = self.library.get_ch3_hole()
        cl = self.library.get_kappa_cl()

        # Heavy chain 1: Target Fab arm (Knob)
        # VH - CH1 - Hinge - CH2 - CH3(Knob)
        heavy_chain_1 = AntibodyChain(
            name="HC1_target_fab_knob",
            sequence=target_vh + ch1 + hinge + ch2 + ch3_knob,
            chain_type="heavy",
            components=["VH(target)", "CH1", "Hinge", "CH2", "CH3(Knob)"],
        )

        # Heavy chain 2: CD3 VHH arm (Hole)
        # VHH - Hinge - CH2 - CH3(Hole)
        heavy_chain_2 = AntibodyChain(
            name="HC2_CD3_vhh_hole",
            sequence=cd3_binder + hinge + ch2 + ch3_hole,
            chain_type="heavy",
            components=["VHH(CD3)", "Hinge", "CH2", "CH3(Hole)"],
        )

        # Light chain: Target only
        # VL - CL
        light_chain = AntibodyChain(
            name="LC_target",
            sequence=target_vl + cl,
            chain_type="light",
            components=["VL(target)", "CL"],
        )

        return BispecificConstruct(
            name=name,
            format_type=self.format_name,
            chains=[heavy_chain_1, heavy_chain_2, light_chain],
            target_1=target_name,
            target_2="CD3",
            notes=(
                "Asymmetric Fab + VHH format. "
                "Target arm is Fab, CD3 arm is single-domain VHH fused to Fc. "
                "Knob-in-hole for HC heterodimerization. "
                "Smallest asymmetric bispecific format."
            ),
        )


def assemble_fab_vhh(
    target_vh: str,
    target_vl: str,
    cd3_vhh: str,
    name: str = "fab_vhh_bispecific",
    target_name: str = "HER2",
    sequence_library: Optional[SequenceLibrary] = None,
) -> BispecificConstruct:
    """Convenience function to assemble Fab + VHH bispecific.

    Args:
        target_vh: VH sequence for tumor target Fab arm.
        target_vl: VL sequence for tumor target Fab arm.
        cd3_vhh: VHH sequence for CD3 arm.
        name: Name for the construct.
        target_name: Name of tumor target.
        sequence_library: Custom sequence library.

    Returns:
        BispecificConstruct with 3 chains.
    """
    formatter = FabVhhFormatter(sequence_library)
    return formatter.assemble(
        target_vh=target_vh,
        target_vl=target_vl,
        cd3_binder=cd3_vhh,
        name=name,
        target_name=target_name,
    )
