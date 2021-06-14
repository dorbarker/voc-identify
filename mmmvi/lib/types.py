import pysam

from typing import Dict, List, Tuple, Optional

Position = Tuple[int, ...]
Mutation = Tuple[Optional[str], ...]
Mutations = Dict[Position, Mutation]
VoCs = Dict[str, Mutations]
Reads = List[pysam.AlignedSegment]
MutationResults = Dict[str, List[Position]]
VoCResults = Dict[str, MutationResults]
