#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - stitch adapter prototype
"""

from enum import Enum
from typing import Optional

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

class STITCHEdgeType(Enum):
    """
    STITCH edge types.
    """
    METABOLITE_TO_PROTEIN = "METABOLITE_TO_PROTEIN"

class STITCHMetaboliteTOProteinEdgeField(Enum):
    """
    STITCH metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"
    REACTION_ID = "REACTION_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    MODE = "MODE"
    DATABASE = "DATABASE"
    EXPERIMENT = "EXPERIMENT"
    PREDICTION = "PREDICTION"
    TEXTMINING = "TEXTMINING"
    COMBINED_SCORE = "COMBINED_SCORE"


class STITCHAdapter:









