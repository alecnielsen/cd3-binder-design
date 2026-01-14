"""CrossMab bispecific format assembly.

CrossMab uses CH1-CL domain swap on one arm to ensure correct light chain pairing.
Both arms are Fab format with knob-in-hole for heavy chain heterodimerization.

Structure:
- Heavy chain 1 (Knob): VH(target) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2 (Hole): VH(CD3) - CL - Hinge - CH2 - CH3(Hole)  [CrossMab swap]
- Light chain 1: VL(target) - CL
- Light chain 2: VL(CD3) - CH1  [CrossMab swap]
"""

from typing import Optional

from src.formatting.base import (
    BispecificFormatter,
    BispecificConstruct,
    AntibodyChain,
    SequenceLibrary,
)


class CrossMabFormatter(BispecificFormatter):
    """Formatter for CrossMab (Fab x Fab) bispecifics."""

    @property
    def format_name(self) -> str:
        return "crossmab"

    def assemble(
        self,
        target_vh: str,
        target_vl: str,
        cd3_binder: str,  # VH for CD3
        cd3_binder_vl: Optional[str] = None,  # VL for CD3 (required for CrossMab)
        name: str = "crossmab_bispecific",
        target_name: str = "HER2",
    ) -> BispecificConstruct:
        """Assemble CrossMab bispecific.

        Args:
            target_vh: VH sequence for tumor target arm.
            target_vl: VL sequence for tumor target arm.
            cd3_binder: VH sequence for CD3 arm.
            cd3_binder_vl: VL sequence for CD3 arm (required).
            name: Name for the construct.
            target_name: Name of tumor target.

        Returns:
            BispecificConstruct with 4 chains.
        """
        if cd3_binder_vl is None:
            raise ValueError("CrossMab requires VL for CD3 arm. Use Fab+VHH format for VHH binders.")

        # Get constant region sequences
        ch1 = self.library.get_ch1()
        hinge = self.library.get_hinge()
        ch2 = self.library.get_ch2()
        ch3_knob = self.library.get_ch3_knob()
        ch3_hole = self.library.get_ch3_hole()
        cl = self.library.get_kappa_cl()

        # Heavy chain 1: Target arm (standard Fab, Knob)
        # VH - CH1 - Hinge - CH2 - CH3(Knob)
        heavy_chain_1 = AntibodyChain(
            name="HC1_target_knob",
            sequence=target_vh + ch1 + hinge + ch2 + ch3_knob,
            chain_type="heavy",
            components=["VH(target)", "CH1", "Hinge", "CH2", "CH3(Knob)"],
        )

        # Heavy chain 2: CD3 arm (CrossMab swap, Hole)
        # VH - CL - Hinge - CH2 - CH3(Hole)
        # Note: CH1 is replaced with CL for CrossMab
        heavy_chain_2 = AntibodyChain(
            name="HC2_CD3_hole_crossmab",
            sequence=cd3_binder + cl + hinge + ch2 + ch3_hole,
            chain_type="heavy",
            components=["VH(CD3)", "CL(CrossMab)", "Hinge", "CH2", "CH3(Hole)"],
        )

        # Light chain 1: Target arm (standard)
        # VL - CL
        light_chain_1 = AntibodyChain(
            name="LC1_target",
            sequence=target_vl + cl,
            chain_type="light",
            components=["VL(target)", "CL"],
        )

        # Light chain 2: CD3 arm (CrossMab swap)
        # VL - CH1
        # Note: CL is replaced with CH1 for CrossMab
        light_chain_2 = AntibodyChain(
            name="LC2_CD3_crossmab",
            sequence=cd3_binder_vl + ch1,
            chain_type="light",
            components=["VL(CD3)", "CH1(CrossMab)"],
        )

        return BispecificConstruct(
            name=name,
            format_type=self.format_name,
            chains=[heavy_chain_1, heavy_chain_2, light_chain_1, light_chain_2],
            target_1=target_name,
            target_2="CD3",
            notes=(
                "CrossMab format with CH1-CL swap on CD3 arm for correct LC pairing. "
                "Knob-in-hole mutations for HC heterodimerization. "
                "LALA mutations in CH2 for Fc silencing."
            ),
        )


def assemble_crossmab(
    target_vh: str,
    target_vl: str,
    cd3_vh: str,
    cd3_vl: str,
    name: str = "crossmab_bispecific",
    target_name: str = "HER2",
    sequence_library: Optional[SequenceLibrary] = None,
) -> BispecificConstruct:
    """Convenience function to assemble CrossMab bispecific.

    Args:
        target_vh: VH sequence for tumor target arm.
        target_vl: VL sequence for tumor target arm.
        cd3_vh: VH sequence for CD3 arm.
        cd3_vl: VL sequence for CD3 arm.
        name: Name for the construct.
        target_name: Name of tumor target.
        sequence_library: Custom sequence library.

    Returns:
        BispecificConstruct with 4 chains.
    """
    formatter = CrossMabFormatter(sequence_library)
    return formatter.assemble(
        target_vh=target_vh,
        target_vl=target_vl,
        cd3_binder=cd3_vh,
        cd3_binder_vl=cd3_vl,
        name=name,
        target_name=target_name,
    )
