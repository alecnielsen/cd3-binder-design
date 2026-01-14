"""IgG-VHH Morrison format bispecific assembly.

Symmetric format with full IgG (bivalent for tumor target) and VHH fused
to C-terminus of both heavy chains (bivalent for CD3).

Structure:
- Heavy chain: VH(target) - CH1 - Hinge - CH2 - CH3 - Linker - VHH(CD3)
- Light chain: VL(target) - CL

Valency:
- 2x tumor target binding (bivalent Fab)
- 2x CD3 binding (bivalent VHH)

Note: This is a SYMMETRIC format - no knob-in-hole needed.
Smallest Morrison format due to single-domain VHH.
"""

from typing import Optional

from src.formatting.base import (
    BispecificFormatter,
    BispecificConstruct,
    AntibodyChain,
    SequenceLibrary,
)


class IggVhhFormatter(BispecificFormatter):
    """Formatter for IgG-(VHH)2 Morrison bispecifics."""

    @property
    def format_name(self) -> str:
        return "igg_vhh"

    def assemble(
        self,
        target_vh: str,
        target_vl: str,
        cd3_binder: str,  # VHH for CD3
        cd3_binder_vl: Optional[str] = None,  # Ignored for VHH
        name: str = "igg_vhh_bispecific",
        target_name: str = "HER2",
        fusion_linker: Optional[str] = None,
    ) -> BispecificConstruct:
        """Assemble IgG-VHH bispecific (Morrison format).

        Args:
            target_vh: VH sequence for tumor target.
            target_vl: VL sequence for tumor target.
            cd3_binder: VHH sequence for CD3.
            cd3_binder_vl: Ignored (VHH has no VL).
            name: Name for the construct.
            target_name: Name of tumor target.
            fusion_linker: Custom Fc-to-VHH fusion linker.

        Returns:
            BispecificConstruct with 2 chains (symmetric HC, one LC type).
        """
        # Get constant region sequences
        ch1 = self.library.get_ch1()
        hinge = self.library.get_hinge()
        ch2 = self.library.get_ch2()
        ch3 = self.library.get_ch3_standard()  # Standard CH3, no knob-hole (symmetric)
        cl = self.library.get_kappa_cl()

        # Fusion linker (shorter for VHH than scFv)
        fusion_link = fusion_linker or self.library.get_fc_fusion_linker()

        # Heavy chain: VH - CH1 - Hinge - CH2 - CH3 - Linker - VHH
        # Both heavy chains are identical (symmetric)
        heavy_chain = AntibodyChain(
            name="HC_target_igg_vhh",
            sequence=target_vh + ch1 + hinge + ch2 + ch3 + fusion_link + cd3_binder,
            chain_type="heavy",
            components=["VH(target)", "CH1", "Hinge", "CH2", "CH3", f"Linker({len(fusion_link)}aa)", "VHH(CD3)"],
        )

        # Light chain: VL - CL
        light_chain = AntibodyChain(
            name="LC_target",
            sequence=target_vl + cl,
            chain_type="light",
            components=["VL(target)", "CL"],
        )

        return BispecificConstruct(
            name=name,
            format_type=self.format_name,
            chains=[heavy_chain, light_chain],
            target_1=target_name,
            target_2="CD3",
            notes=(
                "IgG-(VHH)2 Morrison format (symmetric). "
                "Bivalent for both targets: 2x Fab(target), 2x VHH(CD3). "
                "VHH fused to C-terminus of both heavy chains. "
                f"Fc fusion linker: {fusion_link}"
            ),
        )


def assemble_igg_vhh(
    target_vh: str,
    target_vl: str,
    cd3_vhh: str,
    name: str = "igg_vhh_bispecific",
    target_name: str = "HER2",
    fusion_linker: Optional[str] = None,
    sequence_library: Optional[SequenceLibrary] = None,
) -> BispecificConstruct:
    """Convenience function to assemble IgG-VHH bispecific.

    Args:
        target_vh: VH sequence for tumor target.
        target_vl: VL sequence for tumor target.
        cd3_vhh: VHH sequence for CD3.
        name: Name for the construct.
        target_name: Name of tumor target.
        fusion_linker: Custom Fc-to-VHH fusion linker.
        sequence_library: Custom sequence library.

    Returns:
        BispecificConstruct with 2 chain types.
    """
    formatter = IggVhhFormatter(sequence_library)
    return formatter.assemble(
        target_vh=target_vh,
        target_vl=target_vl,
        cd3_binder=cd3_vhh,
        name=name,
        target_name=target_name,
        fusion_linker=fusion_linker,
    )
