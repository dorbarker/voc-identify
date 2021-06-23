import pysam

from typing import Dict, List, Tuple, Optional, Set, Union

Position = Tuple[int, ...]
Mutation = Tuple[Optional[str], ...]
Mutations = Dict[Position, Mutation]
VoCs = Dict[str, Mutations]
Reads = List[pysam.AlignedSegment]
Reads = Dict[
    str,
    Dict[
        str,
        Union[Set[str], List[Tuple[Optional[int], Optional[int]]], List[Optional[int]]],
    ],
]
MutationResults = Dict[str, List[Position]]
VoCResults = Dict[str, MutationResults]
