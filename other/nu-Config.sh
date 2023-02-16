
# Basic code operation

#    LEAN     # if turned on, all particles of the same type are assumed to be of the same mass.

    PERIODIC
    SELFGRAVITY
    #RANDOMIZE_DOMAINCENTER
    
# Gravity options

    PMGRID=128
    #TREEPM_NOTIMESPLIT
    #ASMTH=0.01
    MULTIPOLE_ORDER=1     # multipole order for the tree force expansion. Stock is 2, best to set it to 1.
#    TREE_NUM_BEFORE_NODESPLIT=6
    
# Softening types and particle types

    NSOFTCLASSES=6            # set to maximum number of type you expect to use
    NTYPES=6

# Floating point accuracy

    POSITIONS_IN_32BIT
    DOUBLEPRECISION=2

# Miscellaneous code options

    POWERSPEC_ON_OUTPUT
#    CB_VELDIV

#    GADGET2_HEADER 
# IC generation via N-GenIC

#    FOF
#    FOF_PRIMARY_LINK_TYPES=2
#    FOF_SECONDARY_LINK_TYPES=4
#    FOF_GROUP_MIN_LEN=32
#    FOF_LINKLENGTH=0.2

#    SUBFIND
#    SUBFIND_STORE_LOCAL_DENSITY

#Hybrid settings
    NGENIC=128
# turn on these if running Hybrid
#    ADDITIONAL_GRID
#    THERMAL_VEL_IC
#    CB_PHASE

#    MFLR_RST

# turn on this if running Linear response
    CREATE_GRID




#    IDS_32BIT
#    NGENIC_TEST
