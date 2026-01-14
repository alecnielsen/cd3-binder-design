"""IgG-scFv Morrison format bispecific assembly.

Symmetric format with full IgG (bivalent for tumor target) and scFv fused
to C-terminus of both heavy chains (bivalent for CD3).

Structure:
- Heavy chain: VH(target) - CH1 - Hinge - CH2 - CH3 - Linker - scFv(CD3)
- Light chain: VL(target) - CL

Valency:
- 2x tumor target binding (bivalent Fab)
- 2x CD3 binding (bivalent scFv)

Note: This is a SYMMETRIC format - no knob-in-hole needed.
"""

from typing import Optional

from src.formatting.base import (
    BispecificFormatter,
    BispecificConstruct,
    AntibodyChain,
    SequenceLibrary,
)


class IggScfvFormatter(BispecificFormatter):
    """Formatter for IgG-(scFv)2 Morrison bispecifics."""

    @property
    def format_name(self) -> str:
        return "igg_scfv"

    def assemble(
        self,
        target_vh: str,
        target_vl: str,
        cd3_binder: str,  # VH for CD3 scFv
        cd3_binder_vl: Optional[str] = None,  # VL for CD3 scFv
        name: str = "igg_scfv_bispecific",
        target_name: str = "HER2",
        scfv_linker: Optional[str] = None,
        fusion_linker: Optional[str] = None,
    ) -> BispecificConstruct:
        """Assemble IgG-scFv bispecific (Morrison format).

        Args:
            target_vh: VH sequence for tumor target.
            target_vl: VL sequence for tumor target.
            cd3_binder: VH sequence for CD3 scFv.
            cd3_binder_vl: VL sequence for CD3 scFv (required).
            name: Name for the construct.
            target_name: Name of tumor target.
            scfv_linker: Custom scFv internal linker.
            fusion_linker: Custom Fc-to-scFv fusion linker.

        Returns:
            BispecificConstruct with 2 chains (symmetric HC, one LC type).
        """
        if cd3_binder_vl is None:
            raise ValueError("IgG-scFv format requires VL for CD3 scFv. Use IgG-VHH format for VHH binders.")

        # Get constant region sequences
        ch1 = self.library.get_ch1()
        hinge = self.library.get_hinge()
        ch2 = self.library.get_ch2()
        ch3 = self.library.get_ch3_standard()  # Standard CH3, no knob-hole (symmetric)
        cl = self.library.get_kappa_cl()

        # Linkers
        scfv_link = scfv_linker or self.library.get_scfv_linker()
        fusion_link = fusion_linker or self.library.get_fc_fusion_linker()

        # Create scFv from CD3 VH/VL
        cd3_scfv = self.make_scfv(cd3_binder, cd3_binder_vl, scfv_link)

        # Heavy chain: VH - CH1 - Hinge - CH2 - CH3 - Linker - scFv
        # Both heavy chains are identical (symmetric)
        heavy_chain = AntibodyChain(
            name="HC_target_igg_scfv",
            sequence=target_vh + ch1 + hinge + ch2 + ch3 + fusion_link + cd3_scfv,
            chain_type="heavy",
            components=["VH(target)", "CH1", "Hinge", "CH2", "CH3", f"Linker({len(fusion_link)}aa)", "scFv(CD3)"],
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
                "IgG-(scFv)2 Morrison format (symmetric). "
                "Bivalent for both targets: 2x Fab(target), 2x scFv(CD3). "
                "scFv fused to C-terminus of both heavy chains. "
                f"scFv linker: {scfv_link}, Fc fusion linker: {fusion_link}"
            ),
        )


def assemble_igg_scfv(
    target_vh: str,
    target_vl: str,
    cd3_vh: str,
    cd3_vl: str,
    name: str = "igg_scfv_bispecific",
    target_name: str = "HER2",
    scfv_linker: Optional[str] = None,
    fusion_linker: Optional[str] = None,
    sequence_library: Optional[SequenceLibrary] = None,
) -> BispecificConstruct:
    """Convenience function to assemble IgG-scFv bispecific.

    Args:
        target_vh: VH sequence for tumor target.
        target_vl: VL sequence for tumor target.
        cd3_vh: VH sequence for CD3 scFv.
        cd3_vl: VL sequence for CD3 scFv.
        name: Name for the construct.
        target_name: Name of tumor target.
        scfv_linker: Custom scFv internal linker.
        fusion_linker: Custom Fc-to-scFv fusion linker.
        sequence_library: Custom sequence library.

    Returns:
        BispecificConstruct with 2 chain types.
    """
    formatter = IggScfvFormatter(sequence_library)
    return formatter.assemble(
        target_vh=target_vh,
        target_vl=target_vl,
        cd3_binder=cd3_vh,
        cd3_binder_vl=cd3_vl,
        name=name,
        target_name=target_name,
        scfv_linker=scfv_linker,
        fusion_linker=fusion_linker,
    )
