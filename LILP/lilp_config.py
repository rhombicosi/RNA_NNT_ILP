from enum import Enum

VALID_PAIRS = ['AU','UA','CG','GC','GU','UG']
MIN_D = 3
MFE = -1000
SCALE = 100
M = 10000

class LoopType(Enum):
    HAIRPIN = "hairpin"
    INTERNAL = "internal"    
    STEM = "stem"
    BULGE = "bulge"
    MULTI = "multi"

class InternalType(Enum):
    INT11 = "int11"
    INT12 = "int12"
    INT21 = "int21"
    INT1N = "int1n"
    INT22 = "int22"
    INT23 = "int23"
    INT32 = "int32"
    INTGEN = "general"

MAX_LOOP_SIZES = {
        LoopType.HAIRPIN: 10,
        LoopType.INTERNAL: 10,
        LoopType.BULGE: 10,
        LoopType.MULTI: 10,
    }

MAX_NUM_OF_LOOPS = {
        LoopType.HAIRPIN: 1,
        LoopType.INTERNAL: 2,
        LoopType.BULGE: 0,
        LoopType.MULTI: 0,
    }
