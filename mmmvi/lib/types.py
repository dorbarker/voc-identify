import pysam

from typing import Dict, List, Tuple, Optional, Set, Union

Position = Tuple[int, ...]
Mutation = Tuple[Optional[str], ...]
Mutations = Dict[Position, Mutation]
VoCs = Dict[str, Mutations]
Reads = Dict[
    str, Dict[str, Union[Set[str], pysam.AlignedSegment]],
]
MutationResults = Dict[str, List[Position]]
VoCResults = Dict[str, MutationResults]
