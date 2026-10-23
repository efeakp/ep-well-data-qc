"""File readers. Each returns data plus a list of :class:`wellqc.model.Issue` found while parsing."""

from .checkshot import read_checkshot
from .logs import read_las
from .petrel import read_deviation, read_grid_extent, read_well_head, read_well_tops
from .segy import scan_segy
from .sniff import sniff

__all__ = ["read_checkshot", "read_las", "read_deviation", "read_grid_extent", "read_well_head",
           "read_well_tops", "scan_segy", "sniff"]
