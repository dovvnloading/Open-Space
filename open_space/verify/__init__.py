"""Deterministic verification of the repository's own claims.

``python -m open_space.verify`` re-checks, from the committed bytes alone, what
the result folders assert: that recorded checksums match, that links and catalog
counts hold, that the numbers printed in each README exist in the deposited
evidence, that declared claims can be recomputed independently, and that the
documented replay commands still reproduce the committed outputs.

Nothing here judges whether a scientific interpretation is right. Each check
answers a narrower question that has exactly one correct answer.
"""

from .model import FAIL, INFO, PASS, WARN, Finding
from .tree import RepoTree, TreeError

__all__ = ["FAIL", "INFO", "PASS", "WARN", "Finding", "RepoTree", "TreeError"]
