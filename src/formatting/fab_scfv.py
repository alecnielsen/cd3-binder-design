"""Asymmetric Fab + scFv bispecific format assembly.

One arm is a standard Fab (for tumor target), the other arm is an scFv
(for CD3) fused directly to the Fc. Uses knob-in-hole for heterodimerization.

Structure:
- Heavy chain 1 (Knob): VH(target) - CH1 - Hinge - CH2 - CH3(Knob)
- Heavy chain 2 (Hole): scFv(CD3) - Hinge - CH2 - CH3(Hole)
- Light chain: VL(target) - CL

Advantages:
- Simpler than CrossMab (only one light chain)
- Smaller than full IgG formats
"""

from typing import Optional

from src.formatting.base import (
    BispecificFormatter,
    BispecificConstruct,
    AntibodyChain,
    SequenceLibrary,
)


class FabScFvFormatter(BispecificFormatter):
    """Formatter for asymmetric Fab + scFv bispecifics."""

    @property
    def format_name(self) -> str:
        return "fab_scfv"

    def assemble(
        self,
        target_vh: str,
        target_vl: str,
        cd3_binder: str,  # VH for CD3 scFv
        cd3_binder_vl: Optional[str] = None,  # VL for CD3 scFv
        name: str = "fab_scfv_bispecific",
        target_name: str = "HER2",
        scfv_linker: Optional[str] = None,
    ) -> BispecificConstruct:
        """Assemble Fab + scFv bispecific.

        Args:
            target_vh: VH sequence for tumor target Fab arm.
            target_vl: VL sequence for tumor target Fab arm.
            cd3_binder: VH sequence for CD3 scFv.
            cd3_binder_vl: VL sequence for CD3 scFv (required).
            name: Name for the construct.
            target_name: Name of tumor target.
            scfv_linker: Custom scFv linker sequence.

        Returns:
            BispecificConstruct with 3 chains.
        """
        if cd3_binder_vl is None:
            raise ValueError("Fab+scFv format requires VL for CD3 scFv. Use Fab+VHH format for VHH binders.")

        # Get constant region sequences
        ch1 = self.library.get_ch1()
        hinge = self.library.get_hinge()
        ch2 = self.library.get_ch2()
        ch3_knob = self.library.get_ch3_knob()
        ch3_hole = self.library.get_ch3_hole()
        cl = self.library.get_kappa_cl()
        linker = scfv_linker or self.library.get_scfv_linker()

        # Create scFv from CD3 VH/VL
        cd3_scfv = self.make_scfv(cd3_binder, cd3_binder_vl, linker)

        # Heavy chain 1: Target Fab arm (Knob)
        # VH - CH1 - Hinge - CH2 - CH3(Knob)
        heavy_chain_1 = AntibodyChain(
            name="HC1_target_fab_knob",
            sequence=target_vh + ch1 + hinge + ch2 + ch3_knob,
            chain_type="heavy",
            components=["VH(target)", "CH1", "Hinge", "CH2", "CH3(Knob)"],
        )

        # Heavy chain 2: CD3 scFv arm (Hole)
        # scFv - Hinge - CH2 - CH3(Hole)
        heavy_chain_2 = AntibodyChain(
            name="HC2_CD3_scfv_hole",
            sequence=cd3_scfv + hinge + ch2 + ch3_hole,
            chain_type="heavy",
            components=["scFv(CD3)", "Hinge", "CH2", "CH3(Hole)"],
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
                "Asymmetric Fab + scFv format. "
                "Target arm is Fab, CD3 arm is scFv fused to Fc. "
                "Knob-in-hole for HC heterodimerization. "
                f"scFv linker: {linker}"
            ),
        )


def assemble_fab_scfv(
    target_vh: str,
    target_vl: str,
    cd3_vh: str,
    cd3_vl: str,
    name: str = "fab_scfv_bispecific",
    target_name: str = "HER2",
    scfv_linker: Optional[str] = None,
    sequence_library: Optional[SequenceLibrary] = None,
) -> BispecificConstruct:
    """Convenience function to assemble Fab + scFv bispecific.

    Args:
        target_vh: VH sequence for tumor target Fab arm.
        target_vl: VL sequence for tumor target Fab arm.
        cd3_vh: VH sequence for CD3 scFv.
        cd3_vl: VL sequence for CD3 scFv.
        name: Name for the construct.
        target_name: Name of tumor target.
        scfv_linker: Custom scFv linker sequence.
        sequence_library: Custom sequence library.

    Returns:
        BispecificConstruct with 3 chains.
    """
    formatter = FabScFvFormatter(sequence_library)
    return formatter.assemble(
        target_vh=target_vh,
        target_vl=target_vl,
        cd3_binder=cd3_vh,
        cd3_binder_vl=cd3_vl,
        name=name,
        target_name=target_name,
        scfv_linker=scfv_linker,
    )
