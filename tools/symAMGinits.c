#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <blas.h>
#include <janus.h>

#include <ilupackmacros.h>

#define MAX_LINE        255
#define STDERR          stdout
#define STDOUT          stdout
#define PRINT_INFO
#define MAX(A,B)        (((A)>(B))?(A):(B))
#define MIN(A,B)        (((A)<(B))?(A):(B))
#define RESTOL_FUNC(A)  sqrt(A)





#if !defined _DOUBLE_REAL_ && !defined _SINGLE_REAL_
#define _COMPLEX_SYMMETRIC_



#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
#define CONJG(A)      (A)

#ifdef _DOUBLE_REAL_

#define PERMSMWMAMF          DSYMperm_mwm_amf
#define PERMSMWMAMD          DSYMperm_mwm_amd
#define PERMSMWMMMD          DSYMperm_mwm_mmd
#define PERMSMWMRCM          DSYMperm_mwm_rcm
#define PERMSMWMMETISE       DSYMperm_mwm_metis_e
#define PERMSMWMMETISN       DSYMperm_mwm_metis_n

#define PERMSMC64AMF          DSYMperm_mc64_amf
#define PERMSMC64AMD          DSYMperm_mc64_amd
#define PERMSMC64MMD          DSYMperm_mc64_mmd
#define PERMSMC64RCM          DSYMperm_mc64_rcm
#define PERMSMC64METISE       DSYMperm_mc64_metis_e
#define PERMSMC64METISN       DSYMperm_mc64_metis_n

#define PERMSMATCHINGAMF          DSYMperm_matching_amf
#define PERMSMATCHINGAMD          DSYMperm_matching_amd
#define PERMSMATCHINGMMD          DSYMperm_matching_mmd
#define PERMSMATCHINGRCM          DSYMperm_matching_rcm
#define PERMSMATCHINGMETISE       DSYMperm_matching_metis_e
#define PERMSMATCHINGMETISN       DSYMperm_matching_metis_n

#define PERMSMWMAMFFCV        DSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        DSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        DSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        DSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     DSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     DSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       DSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       DSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       DSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       DSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    DSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    DSYMperm_mc64_metis_n_fcv

#define PERMSMATCHINGAMFFCV       DSYMperm_matching_amf_fcv
#define PERMSMATCHINGAMDFCV       DSYMperm_matching_amd_fcv
#define PERMSMATCHINGMMDFCV       DSYMperm_matching_mmd_fcv
#define PERMSMATCHINGRCMFCV       DSYMperm_matching_rcm_fcv
#define PERMSMATCHINGMETISEFCV    DSYMperm_matching_metis_e_fcv
#define PERMSMATCHINGMETISNFCV    DSYMperm_matching_metis_n_fcv

#define PERMSMWMAMFFC        DSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        DSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        DSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        DSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     DSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     DSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       DSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       DSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       DSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       DSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    DSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    DSYMperm_mc64_metis_n_fc

#define PERMSMATCHINGAMFFC       DSYMperm_matching_amf_fc
#define PERMSMATCHINGAMDFC       DSYMperm_matching_amd_fc
#define PERMSMATCHINGMMDFC       DSYMperm_matching_mmd_fc
#define PERMSMATCHINGRCMFC       DSYMperm_matching_rcm_fc
#define PERMSMATCHINGMETISEFC    DSYMperm_matching_metis_e_fc
#define PERMSMATCHINGMETISNFC    DSYMperm_matching_metis_n_fc

#else

#define PERMSMWMAMF          SSYMperm_mwm_amf
#define PERMSMWMAMD          SSYMperm_mwm_amd
#define PERMSMWMMMD          SSYMperm_mwm_mmd
#define PERMSMWMRCM          SSYMperm_mwm_rcm
#define PERMSMWMMETISE       SSYMperm_mwm_metis_e
#define PERMSMWMMETISN       SSYMperm_mwm_metis_n

#define PERMSMC64AMF          SSYMperm_mc64_amf
#define PERMSMC64AMD          SSYMperm_mc64_amd
#define PERMSMC64MMD          SSYMperm_mc64_mmd
#define PERMSMC64RCM          SSYMperm_mc64_rcm
#define PERMSMC64METISE       SSYMperm_mc64_metis_e
#define PERMSMC64METISN       SSYMperm_mc64_metis_n

#define PERMSMATCHINGAMF          SSYMperm_matching_amf
#define PERMSMATCHINGAMD          SSYMperm_matching_amd
#define PERMSMATCHINGMMD          SSYMperm_matching_mmd
#define PERMSMATCHINGRCM          SSYMperm_matching_rcm
#define PERMSMATCHINGMETISE       SSYMperm_matching_metis_e
#define PERMSMATCHINGMETISN       SSYMperm_matching_metis_n

#define PERMSMWMAMFFCV        SSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        SSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        SSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        SSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     SSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     SSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       SSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       SSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       SSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       SSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    SSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    SSYMperm_mc64_metis_n_fcv

#define PERMSMATCHINGAMFFCV       SSYMperm_matching_amf_fcv
#define PERMSMATCHINGAMDFCV       SSYMperm_matching_amd_fcv
#define PERMSMATCHINGMMDFCV       SSYMperm_matching_mmd_fcv
#define PERMSMATCHINGRCMFCV       SSYMperm_matching_rcm_fcv
#define PERMSMATCHINGMETISEFCV    SSYMperm_matching_metis_e_fcv
#define PERMSMATCHINGMETISNFCV    SSYMperm_matching_metis_n_fcv

#define PERMSMWMAMFFC        SSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        SSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        SSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        SSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     SSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     SSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       SSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       SSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       SSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       SSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    SSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    SSYMperm_mc64_metis_n_fc

#define PERMSMATCHINGAMFFC       SSYMperm_matching_amf_fc
#define PERMSMATCHINGAMDFC       SSYMperm_matching_amd_fc
#define PERMSMATCHINGMMDFC       SSYMperm_matching_mmd_fc
#define PERMSMATCHINGRCMFC       SSYMperm_matching_rcm_fc
#define PERMSMATCHINGMETISEFC    SSYMperm_matching_metis_e_fc
#define PERMSMATCHINGMETISNFC    SSYMperm_matching_metis_n_fc

#endif




#else



#ifdef _COMPLEX_SYMMETRIC_
#define CONJG(A)     (A)

#ifdef _SINGLE_COMPLEX_
#define PERMSMWMAMF          CSYMperm_mwm_amf
#define PERMSMWMAMD          CSYMperm_mwm_amd
#define PERMSMWMMMD          CSYMperm_mwm_mmd
#define PERMSMWMRCM          CSYMperm_mwm_rcm
#define PERMSMWMMETISE       CSYMperm_mwm_metis_e
#define PERMSMWMMETISN       CSYMperm_mwm_metis_n

#define PERMSMC64AMF          CSYMperm_mc64_amf
#define PERMSMC64AMD          CSYMperm_mc64_amd
#define PERMSMC64MMD          CSYMperm_mc64_mmd
#define PERMSMC64RCM          CSYMperm_mc64_rcm
#define PERMSMC64METISE       CSYMperm_mc64_metis_e
#define PERMSMC64METISN       CSYMperm_mc64_metis_n

#define PERMSMATCHINGAMF          CSYMperm_matching_amf
#define PERMSMATCHINGAMD          CSYMperm_matching_amd
#define PERMSMATCHINGMMD          CSYMperm_matching_mmd
#define PERMSMATCHINGRCM          CSYMperm_matching_rcm
#define PERMSMATCHINGMETISE       CSYMperm_matching_metis_e
#define PERMSMATCHINGMETISN       CSYMperm_matching_metis_n

#define PERMSMWMAMFFCV        CSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        CSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        CSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        CSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     CSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     CSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       CSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       CSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       CSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       CSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    CSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    CSYMperm_mc64_metis_n_fcv

#define PERMSMATCHINGAMFFCV       CSYMperm_matching_amf_fcv
#define PERMSMATCHINGAMDFCV       CSYMperm_matching_amd_fcv
#define PERMSMATCHINGMMDFCV       CSYMperm_matching_mmd_fcv
#define PERMSMATCHINGRCMFCV       CSYMperm_matching_rcm_fcv
#define PERMSMATCHINGMETISEFCV    CSYMperm_matching_metis_e_fcv
#define PERMSMATCHINGMETISNFCV    CSYMperm_matching_metis_n_fcv

#define PERMSMWMAMFFC        CSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        CSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        CSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        CSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     CSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     CSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       CSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       CSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       CSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       CSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    CSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    CSYMperm_mc64_metis_n_fc

#define PERMSMATCHINGAMFFC       CSYMperm_matching_amf_fc
#define PERMSMATCHINGAMDFC       CSYMperm_matching_amd_fc
#define PERMSMATCHINGMMDFC       CSYMperm_matching_mmd_fc
#define PERMSMATCHINGRCMFC       CSYMperm_matching_rcm_fc
#define PERMSMATCHINGMETISEFC    CSYMperm_matching_metis_e_fc
#define PERMSMATCHINGMETISNFC    CSYMperm_matching_metis_n_fc


#else
#define PERMSMWMAMF          ZSYMperm_mwm_amf
#define PERMSMWMAMD          ZSYMperm_mwm_amd
#define PERMSMWMMMD          ZSYMperm_mwm_mmd
#define PERMSMWMRCM          ZSYMperm_mwm_rcm
#define PERMSMWMMETISE       ZSYMperm_mwm_metis_e
#define PERMSMWMMETISN       ZSYMperm_mwm_metis_n

#define PERMSMC64AMF          ZSYMperm_mc64_amf
#define PERMSMC64AMD          ZSYMperm_mc64_amd
#define PERMSMC64MMD          ZSYMperm_mc64_mmd
#define PERMSMC64RCM          ZSYMperm_mc64_rcm
#define PERMSMC64METISE       ZSYMperm_mc64_metis_e
#define PERMSMC64METISN       ZSYMperm_mc64_metis_n

#define PERMSMATCHINGAMF          ZSYMperm_matching_amf
#define PERMSMATCHINGAMD          ZSYMperm_matching_amd
#define PERMSMATCHINGMMD          ZSYMperm_matching_mmd
#define PERMSMATCHINGRCM          ZSYMperm_matching_rcm
#define PERMSMATCHINGMETISE       ZSYMperm_matching_metis_e
#define PERMSMATCHINGMETISN       ZSYMperm_matching_metis_n

#define PERMSMWMAMFFCV        ZSYMperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        ZSYMperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        ZSYMperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        ZSYMperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     ZSYMperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     ZSYMperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       ZSYMperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       ZSYMperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       ZSYMperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       ZSYMperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    ZSYMperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    ZSYMperm_mc64_metis_n_fcv

#define PERMSMATCHINGAMFFCV       ZSYMperm_matching_amf_fcv
#define PERMSMATCHINGAMDFCV       ZSYMperm_matching_amd_fcv
#define PERMSMATCHINGMMDFCV       ZSYMperm_matching_mmd_fcv
#define PERMSMATCHINGRCMFCV       ZSYMperm_matching_rcm_fcv
#define PERMSMATCHINGMETISEFCV    ZSYMperm_matching_metis_e_fcv
#define PERMSMATCHINGMETISNFCV    ZSYMperm_matching_metis_n_fcv

#define PERMSMWMAMFFC        ZSYMperm_mwm_amf_fc
#define PERMSMWMAMDFC        ZSYMperm_mwm_amd_fc
#define PERMSMWMMMDFC        ZSYMperm_mwm_mmd_fc
#define PERMSMWMRCMFC        ZSYMperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     ZSYMperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     ZSYMperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       ZSYMperm_mc64_amf_fc
#define PERMSMC64AMDFC       ZSYMperm_mc64_amd_fc
#define PERMSMC64MMDFC       ZSYMperm_mc64_mmd_fc
#define PERMSMC64RCMFC       ZSYMperm_mc64_rcm_fc
#define PERMSMC64METISEFC    ZSYMperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    ZSYMperm_mc64_metis_n_fc

#define PERMSMATCHINGAMFFC       ZSYMperm_matching_amf_fc
#define PERMSMATCHINGAMDFC       ZSYMperm_matching_amd_fc
#define PERMSMATCHINGMMDFC       ZSYMperm_matching_mmd_fc
#define PERMSMATCHINGRCMFC       ZSYMperm_matching_rcm_fc
#define PERMSMATCHINGMETISEFC    ZSYMperm_matching_metis_e_fc
#define PERMSMATCHINGMETISNFC    ZSYMperm_matching_metis_n_fc

#endif

#else
#define CONJG(A)             (-(A))

#ifdef _SINGLE_COMPLEX_
#define PERMSMWMAMF          CHERperm_mwm_amf
#define PERMSMWMAMD          CHERperm_mwm_amd
#define PERMSMWMMMD          CHERperm_mwm_mmd
#define PERMSMWMRCM          CHERperm_mwm_rcm
#define PERMSMWMMETISE       CHERperm_mwm_metis_e
#define PERMSMWMMETISN       CHERperm_mwm_metis_n

#define PERMSMC64AMF          CHERperm_mc64_amf
#define PERMSMC64AMD          CHERperm_mc64_amd
#define PERMSMC64MMD          CHERperm_mc64_mmd
#define PERMSMC64RCM          CHERperm_mc64_rcm
#define PERMSMC64METISE       CHERperm_mc64_metis_e
#define PERMSMC64METISN       CHERperm_mc64_metis_n

#define PERMSMATCHINGAMF          CHERperm_matching_amf
#define PERMSMATCHINGAMD          CHERperm_matching_amd
#define PERMSMATCHINGMMD          CHERperm_matching_mmd
#define PERMSMATCHINGRCM          CHERperm_matching_rcm
#define PERMSMATCHINGMETISE       CHERperm_matching_metis_e
#define PERMSMATCHINGMETISN       CHERperm_matching_metis_n

#define PERMSMWMAMFFCV        CHERperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        CHERperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        CHERperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        CHERperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     CHERperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     CHERperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       CHERperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       CHERperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       CHERperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       CHERperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    CHERperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    CHERperm_mc64_metis_n_fcv

#define PERMSMATCHINGAMFFCV       CHERperm_matching_amf_fcv
#define PERMSMATCHINGAMDFCV       CHERperm_matching_amd_fcv
#define PERMSMATCHINGMMDFCV       CHERperm_matching_mmd_fcv
#define PERMSMATCHINGRCMFCV       CHERperm_matching_rcm_fcv
#define PERMSMATCHINGMETISEFCV    CHERperm_matching_metis_e_fcv
#define PERMSMATCHINGMETISNFCV    CHERperm_matching_metis_n_fcv

#define PERMSMWMAMFFC        CHERperm_mwm_amf_fc
#define PERMSMWMAMDFC        CHERperm_mwm_amd_fc
#define PERMSMWMMMDFC        CHERperm_mwm_mmd_fc
#define PERMSMWMRCMFC        CHERperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     CHERperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     CHERperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       CHERperm_mc64_amf_fc
#define PERMSMC64AMDFC       CHERperm_mc64_amd_fc
#define PERMSMC64MMDFC       CHERperm_mc64_mmd_fc
#define PERMSMC64RCMFC       CHERperm_mc64_rcm_fc
#define PERMSMC64METISEFC    CHERperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    CHERperm_mc64_metis_n_fc

#define PERMSMATCHINGAMFFC       CHERperm_matching_amf_fc
#define PERMSMATCHINGAMDFC       CHERperm_matching_amd_fc
#define PERMSMATCHINGMMDFC       CHERperm_matching_mmd_fc
#define PERMSMATCHINGRCMFC       CHERperm_matching_rcm_fc
#define PERMSMATCHINGMETISEFC    CHERperm_matching_metis_e_fc
#define PERMSMATCHINGMETISNFC    CHERperm_matching_metis_n_fc


#else
#define PERMSMWMAMF          ZHERperm_mwm_amf
#define PERMSMWMAMD          ZHERperm_mwm_amd
#define PERMSMWMMMD          ZHERperm_mwm_mmd
#define PERMSMWMRCM          ZHERperm_mwm_rcm
#define PERMSMWMMETISE       ZHERperm_mwm_metis_e
#define PERMSMWMMETISN       ZHERperm_mwm_metis_n

#define PERMSMC64AMF          ZHERperm_mc64_amf
#define PERMSMC64AMD          ZHERperm_mc64_amd
#define PERMSMC64MMD          ZHERperm_mc64_mmd
#define PERMSMC64RCM          ZHERperm_mc64_rcm
#define PERMSMC64METISE       ZHERperm_mc64_metis_e
#define PERMSMC64METISN       ZHERperm_mc64_metis_n

#define PERMSMATCHINGAMF          ZHERperm_matching_amf
#define PERMSMATCHINGAMD          ZHERperm_matching_amd
#define PERMSMATCHINGMMD          ZHERperm_matching_mmd
#define PERMSMATCHINGRCM          ZHERperm_matching_rcm
#define PERMSMATCHINGMETISE       ZHERperm_matching_metis_e
#define PERMSMATCHINGMETISN       ZHERperm_matching_metis_n

#define PERMSMWMAMFFCV        ZHERperm_mwm_amf_fcv
#define PERMSMWMAMDFCV        ZHERperm_mwm_amd_fcv
#define PERMSMWMMMDFCV        ZHERperm_mwm_mmd_fcv
#define PERMSMWMRCMFCV        ZHERperm_mwm_rcm_fcv
#define PERMSMWMMETISEFCV     ZHERperm_mwm_metis_e_fcv
#define PERMSMWMMETISNFCV     ZHERperm_mwm_metis_n_fcv

#define PERMSMC64AMFFCV       ZHERperm_mc64_amf_fcv
#define PERMSMC64AMDFCV       ZHERperm_mc64_amd_fcv
#define PERMSMC64MMDFCV       ZHERperm_mc64_mmd_fcv
#define PERMSMC64RCMFCV       ZHERperm_mc64_rcm_fcv
#define PERMSMC64METISEFCV    ZHERperm_mc64_metis_e_fcv
#define PERMSMC64METISNFCV    ZHERperm_mc64_metis_n_fcv

#define PERMSMATCHINGAMFFCV       ZHERperm_matching_amf_fcv
#define PERMSMATCHINGAMDFCV       ZHERperm_matching_amd_fcv
#define PERMSMATCHINGMMDFCV       ZHERperm_matching_mmd_fcv
#define PERMSMATCHINGRCMFCV       ZHERperm_matching_rcm_fcv
#define PERMSMATCHINGMETISEFCV    ZHERperm_matching_metis_e_fcv
#define PERMSMATCHINGMETISNFCV    ZHERperm_matching_metis_n_fcv

#define PERMSMWMAMFFC        ZHERperm_mwm_amf_fc
#define PERMSMWMAMDFC        ZHERperm_mwm_amd_fc
#define PERMSMWMMMDFC        ZHERperm_mwm_mmd_fc
#define PERMSMWMRCMFC        ZHERperm_mwm_rcm_fc
#define PERMSMWMMETISEFC     ZHERperm_mwm_metis_e_fc
#define PERMSMWMMETISNFC     ZHERperm_mwm_metis_n_fc

#define PERMSMC64AMFFC       ZHERperm_mc64_amf_fc
#define PERMSMC64AMDFC       ZHERperm_mc64_amd_fc
#define PERMSMC64MMDFC       ZHERperm_mc64_mmd_fc
#define PERMSMC64RCMFC       ZHERperm_mc64_rcm_fc
#define PERMSMC64METISEFC    ZHERperm_mc64_metis_e_fc
#define PERMSMC64METISNFC    ZHERperm_mc64_metis_n_fc

#define PERMSMATCHINGAMFFC       ZHERperm_matching_amf_fc
#define PERMSMATCHINGAMDFC       ZHERperm_matching_amd_fc
#define PERMSMATCHINGMMDFC       ZHERperm_matching_mmd_fc
#define PERMSMATCHINGRCMFC       ZHERperm_matching_rcm_fc
#define PERMSMATCHINGMETISEFC    ZHERperm_matching_metis_e_fc
#define PERMSMATCHINGMETISNFC    ZHERperm_matching_metis_n_fc

#endif

#endif



#endif



void SYMAMGINITS(CSRMAT *A, ILUPACKPARAM *param)
{
  /*
    A           complex sammetric sparse matrix in compressed row storage

    param       data structure containing several parameters
                param.ipar[40] miscellaneous integer parameters
                param.fpar[40] miscellaneous real  parameters
                               
    ipar        integer parameters ipar[0,...,33]
                ipar[0]     elbow space factor. Predict the fill-in required
                            by the multilevel ILU by
                            ipar[0]*nnz(A)
                            default: 5
                ipar[1]     
                            
                            
                            
                            
		            
                ipar[2]     
                            
                            
                            
                            
		            
                ipar[3]     lfil parameter. Limit the number of nonzeros per
                            column in L / row in U by ipar[3]
                            default: n+1
                ipar[4]     flag that indicates different types of transposed
                            systems
                            Bit 0:  A^T      instead of A is stored  
                            Bit 1:  CONJG(A) instead of A is stored

                            Bit 2:  solve A^Tx=b      instead of Ax=b  
                            Bit 3:  solve CONJG(A)x=b instead of Ax=b  

                            REMARK: Note that P is derived from A. So if A^T is
			    stored instead of A, it follows that P^T is stored
			    instead of P. Similar relations hold for CONJG(A)

                ipar[5]     reserved for choice of iterative solver 
		            currently unused
			    default: 1

		ipar[6]     flags for the configuration of the multilevel ILU
                            available flags
			    inverse based dropping:                   DROP_INVERSE
			    don't shift away zero pivots(iluc):       NO_SHIFT
			    Tismentsky update:                        TISMENETSKY_SC
			    repeated ILU(iluc):                       REPEAT_FACT
			    improved estimate ||L^{-1}||,||U^{-1}||:  IMPROVED_ESTIMATE
			    diagonal compensation:                    DIAGONAL_COMPENSATION
			    reduce the partial LU to non-coarse
			    part(piluc,mpiluc):                       COARSE_REDUCE
			    use different pivoting strategy, if 
			    the regular reordering fails:             FINAL_PIVOTING
			    initial system should be reordered
			    using an initial strategy:                PREPROCESS_INITIAL_SYSTEM
			    subsequent systems should be reordered
			    using the regular strategy:               PREPROCESS_SUBSYSTEMS
			    use mpiluc as template instead of piluc:  MULTI_PILUC
			    
			    default:  DROP_INVERSE
                                     |PREPROCESS_INITIAL_SYSTEM
				     |PREPROCESS_SUBSYSTEMS
				     |IMPROVED_ESTIMATE
				     |FINAL_PIVOTING
				     |SIMPLE_SC
				     |DISCARD_MATRIX

		ipar[7]     decide which kind of scaling should be combined
                            with the three permutations
			    Bit 0-2 initial preprocessing
                            Bit 0
                                     (no)  left scaling        0/  1
                            Bit 1
                                     (no)  right scaling       0/  2
                            Bit 2    row/column scaling first  0/  4
			    default: 3

			    Bit 3-5 regular reordering
                            Bit 3
                                     (no)  left scaling        0/  8
                            Bit 4
                                     (no)  right scaling       0/ 16
                            Bit 5    row/column scaling first  0/ 32
			    default: 24

			    Bit 6-8 final pivoting
                            Bit 6
                                     (no)  left scaling        0/ 64
                            Bit 7
                                     (no)  right scaling       0/128
                            Bit 8    row/column scaling first  0/256
			    default: 192

                            Bit 9-10 indicate which ordering is currently
                                     in use
                             512     initial preprocessing
                            1024     regular reordering
                            1536     final pivoting
			    default: 0 (no preference)
			    sum default: 3+24+192


		ipar[8]     decide which kind of norm should be used
                            with the three permutations
			    Bit  0- 4 initial preprocessing
                            0         maximum norm
                            1         1-norm
                            2         2-norm 
                            3         square root of the diagonal entries
                            4         symmetric scaling that keeps the entries
			              below one in absolute value
			    default: 4

			    Bit  5- 9 regular reordering
                            0         maximum norm
                            32*1      1-norm
                            32*2      2-norm 
                            32*3      square root of the diagonal entries
                            32*4      symmetric scaling that keeps the entries
			              below one in absolute value
			    default: 32*4

			    Bit 10-14 final pivoting
                            0         maximum norm
                            1024*1    1-norm
                            1024*2    2-norm 
                            1024*3    square root of the diagonal entries
                            1024*4    symmetric scaling that keeps the entries
			              below one in absolute value
			    default: 1024*4
			    sum default: 4+32*4+1024*4

                ipar[9]     lfil parameter for the approximate Schur 
		            complement. Limit the number of nonzeros per
                            row in S by ipar[9]
                            default: n+1

		ipar[10]    define type of AMG
                            0    multilevel ILU (default)
                            1    simple AMG by replacing the solution
                                 of the approximate Schur complement 
                                 recursively by the block ILU from the
				 next level (similar to AMLI)
                            2    traditional AMG 
		ipar[11]    number recursive AMG calls (only if ipar[10]>0)
                            default: 1  (V-cycle)
                                     2  (W-cycle)
		ipar[12]    number of pre-smoothing steps
		            default: 1
		ipar[13]    number of post-smoothing steps
                            default: 1
                ipar[14]    type of pre-smoother
                            default: 1  (block) Gauss-Seidel forward
                                     2  (block) Gauss-Seidel backward
                                     3  (damped block) Jacobi
                                     4  ILU
                ipar[15]    type of post-smoother
                                     1  (block) Gauss-Seidel forward
                            default: 2  (block) Gauss-Seidel backward
                                     3  (damped block) Jacobi
                                     4  ILU

		ipar[16,...,19] currently unused

                ipar[20,...33] parameters used by SPARSKIT solvers
	        ipar[20]    needed to init the iterative solver
		            default: 0
		ipar[21]    preconditioning side (left 1, right 2, both 3)
		            default: 2
		ipar[22]    stopping criteria 
		                     1 (relative residual)
		                     2 (relative ||b||)
			    default: 3 (relative error in the energy norm)
		ipar[23]    workspace size requested by SPARSKIT solver
		            dependent on ipar[5]. Currently only GMRES is
                            supported.
			    default: (n+3)*(ipar[4]+2)+(ipar[4]+1)*ipar[4]/2;

  	        ipar[24]    restart length for GMRES,FOM,... 
		            default: 30
		ipar[25]    maximum number of iteration steps 
		            default: MIN(1.1*n+10,500)
		ipar[26,...33] further parameters used by SPARSKIT solvers
		            default: 0


    fpar        floating point parameters 
                fpar[0]     drop tolerance for triangular factors
                            default:  1e-2
                fpar[1]     drop tolerance for the approximate Schur complement
                            default:  1e-2
                fpar[2]     norm bound for the inverse triangular factors
                            default:  5
                fpar[3]     define when a matrix is dense enough to switch
                            to full matrix processing 
                            fpar[3]*nnz(S) > n^2 
                            default:  3.0
                fpar[4]     Define the minimum block size which allows the
                            the multilevel ILU to continue with the default 
                            reordering strategy (including scaling and 
                            pivoting) before switching to a different strategy
                            if a matrix A of size "n" should be factored at the
                            current level, then piluc yiels up to a permutation
                                 / B  F \
                            A -> |      | such that B is factored
                                 \ E  C /
                            We continue the multilevel factorization if
                            n-nB<=n*fpar[4]
			    default: 0.75
                fpar[5]     Provided that a final static pivoting strategy
                            is requested, we may encounter the situation when
                            the regular reodering strategy becomes ineffective,
                            i.e.
                                 / B  F \
                            A -> |      | such that B is factored
                                 \ E  C /
                            with B of tiny size (n-nB > n*fpar[4]). In this
                            situation we may switch to a different final 
                            pivoting strategy. Like the regular reordering
                            strategy, we may also encounter a situation when
                            the final reordering strategy fails due to
                            a small block B
                            We continue the multilevel factorization with the
                            final pivoting strategy as long as
                            n-nB<=n*fpar[5]
			    default: 0.75
                fpar[6]     ONLY used if MULTI_PILUC is requested. Multiple,
                            repeated use of piluc without switching to a new 
                            level is continued as long as
                            "number of entries skipped" < 
                            "number of entries skipped previously" * fpar[6]
			    default: 0.75

		fpar[7]     drop tolerance for ILDLC, used relative to 
                            fpar[0].	    
			    default: fpar[0]

                fpar[8]     drop tolerance for the remaining Schur complement
                            default: 1e-4 (squared drop tolerance) 
			    	    
                fpar[9]     Define the minimum block size which allows
                            the multilevel ILU to continue with the 
			    current condest value before increasing
                            condest. In principle this is given by
                            fpar[4], except that the condest parameter
                            must be increased in this case. Here the
                            condest parameter is increased to prevent
                            a small leading block on the next level
                            even if the leading block was still large 
                            enough.
			    default: fpar[4]-0.6*(1-fpar[4])

                fpar[10]    for partial factorization define when the number
                            of constraints starts dominating the whole system 

                fpar[11,...,19] currently unused

                fpar[20,...,30] floating parameters required for SPARSKIT 
		            solvers
		fpar[20]    relative error tolerance
		            default: sqrt(MACHEPS)
		fpar[21]    absolute error tolerance
		            default: 0
			    fpar[22]=
		ipar[23,...30] further parameters used by SPARSKIT solvers
		            default: 0

    Authors:

	Matthias Bollhoefer, TU Braunschweig

    Date:

        October 2003, April 2009. ILUPACK V2.3

    Notice:

	Copyright (c) 2009 by TU Braunschweig.  All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Availability:

	This file is located at

	http://ilupack.tu-bs.de/
   */

   float  systime, time_start, time_stop;
   integer i;

   double ILUPACK_secnds_loc[ILUPACK_secnds_length];
   size_t ILUPACK_mem_loc[ILUPACK_mem_length];


   evaluate_time(&time_start,&systime);
   for (i=0; i<ILUPACK_secnds_length; i++) 
       ILUPACK_secnds_loc[i]=0;
   ILUPACK_secnds_loc[8]=time_start;

   // set integer parameters for ILU
   param->ipar[0]=10;
   param->ipar[3]=A->nr+1;
   param->ipar[4]=0; 

   param->ipar[5]=4; // complex symmetric QMR
   
   param->ipar[6]=DROP_INVERSE
          |PREPROCESS_INITIAL_SYSTEM
          |PREPROCESS_SUBSYSTEMS
          |IMPROVED_ESTIMATE
          |FINAL_PIVOTING
          |SIMPLE_SC
          |DISCARD_MATRIX
          |COARSE_REDUCE
          |DECOUPLE_CONSTRAINTSHH;

   param->ipar[7]=3+24+192;
   param->ipar[8]=(1+32+1024)*4;
   param->ipar[9]=A->nr+1;


   // define type of AMG, default: 0  multilevel ILU 
   param->ipar[10]=0;
   // number recursive AMG calls, default: 1  V-cycle
   param->ipar[11]=1;
   // number of pre-smoothing steps
   param->ipar[12]=1;
   // number of post-smoothing steps
   param->ipar[13]=1;
   // type of pre-smoother, default: 1 (block) Gauss-Seidel forward
   param->ipar[14]=1;
   // type of pre-smoother, default: 2 (block) Gauss-Seidel backward
   param->ipar[15]=2;

   // buffer to store 'true' ipar[6] from the time of the factorization
   param->ipar[16]=0;

   // set integer parameters for SPARSKIT solvers
   // init solver
   param->ipar[20]=0;
   // so far use right preconditioning
   param->ipar[21]=2;
   // stopping criteria, backward error
   param->ipar[22]=3;
   // number of restart length for GMRES,FOM,... UNUSED
   param->ipar[24]=30;
   // maximum number of iteration steps 
   param->ipar[25]=1.1*A->nr+10;
   param->ipar[25]=MIN(param->ipar[25],500);
   // init with zeros
   param->niter=0;
   param->ipar[26]=
     param->ipar[27]=
     param->ipar[28]=
     param->ipar[29]=
     param->ipar[30]=
     param->ipar[31]=
     param->ipar[32]=
     param->ipar[33]=0;


   // set floating parameters for ILU
   param->fpar[0]=1.0/10.0; // 1.0;
   param->fpar[1]=1e-2;
   param->fpar[2]=5;
   param->fpar[3]=3.0;
   param->fpar[4]=0.75;
   param->fpar[5]=0.75;
   param->fpar[6]=0.75;
   param->fpar[7]=param->fpar[1]/param->fpar[2];
   param->fpar[8]=1e-3;
   param->fpar[9]=param->fpar[4]-0.6*(1-param->fpar[4]);


   // threshold for partial factorization
   // denote when the number of constraints dominates the whole system
   // default: two third of all unknowns are constraints 
   param->fpar[10]=0.67;

   // adapt the drop tolerance for the remaining Schur complement
   // with respect to the order of the approximate Schur complement
   param->fpar[8]=0.1*param->fpar[1];


   // set floating parameters for SPARSKIT solvers
   param->fpar[20]=sqrt(dgeteps());
#if defined _SINGLE_REAL_ || defined _SINGLE_COMPLEX_   
   /* for single precision sqrt(eps) might by too small to
      be reached. Usually the accuracy for single precision 
      has only half as much digits as double precision, which
      is the square root of eps. To raise the threshold to a 
      reasonable size we take the geometric mean between
      eps^{1/2} and eps^{1/4}.
   */
   param->fpar[20]=sqrt(param->fpar[20]*sqrt(param->fpar[20]));
#endif
   // absolute error tolerance 
   param->fpar[21]=0;
   // init with zeros 
   param->fpar[22]=
     param->fpar[23]=
     param->fpar[24]=
     param->fpar[25]=
     param->fpar[26]=
     param->fpar[27]=
     param->fpar[28]=
     param->fpar[29]=
     param->fpar[30]=0;
       

   param->nibuff=0;
   param->ndbuff=0;
   param->ibuff=NULL;
   param->dbuff=NULL;

   param->nju=0;
   param->njlu=0;
   param->nalu=0;
   param->ju =NULL;
   param->jlu=NULL;
   param->alu=NULL;
   param->ntestvector=0;
   param->testvector=NULL;

   param->niaux=0;
   param->ndaux=0;
   param->iaux=NULL;
   param->daux=NULL;

   param->rcomflag=0;
   param->returnlabel=0;

   // ******************************************************************
   // new parameter selection scheme
   // ******************************************************************
   // turn on matching
   param->matching      =1;
   // preselect permutation strategy
   param->ordering      ="amd";
   // pre-select associated permutation drivers
#ifdef _MC64_MATCHING_
   param->perm0         =PERMSMC64AMD;
   param->perm          =PERMSMC64AMD;
#elif defined _PARDISO_MATCHING_
   param->perm0         =PERMSMWMAMD;
   param->perm          =PERMSMWMAMD;
#else /* _MUMPS_MATCHING_ assumed by default */
   param->perm0         =PERMSMATCHINGAMD;
   param->perm          =PERMSMATCHINGAMD;
#endif
   param->permf         =SPDPERMPP;
   param->ipar[6]|=PREPROCESS_INITIAL_SYSTEM
                  |PREPROCESS_SUBSYSTEMS;
   param->ipar[6]&=~FINAL_PIVOTING;

   param->droptol       =param->fpar[1];
   param->droptolS      =param->fpar[8];
   param->droptolc      =1e4*dgeteps();
   param->condest       =param->fpar[2];
   param->restol        =param->fpar[20];
   param->maxit         =param->ipar[25];
   param->elbow         =(double)param->ipar[0];
   param->lfil          =param->ipar[3];
   param->lfilS         =param->ipar[9];

   // test vector
   if ((param->ipar[6]&DIAGONAL_COMPENSATION) && 
       (param->ipar[6]&STATIC_TESTVECTOR))
      param->typetv        ="static";
   else if ((param->ipar[6]&DIAGONAL_COMPENSATION) && 
	    (param->ipar[6]&DYNAMIC_TESTVECTOR))
      param->typetv        ="dynamic";
   else 
      param->typetv        ="none";
   param->tv            =NULL;
   param->ind           =NULL;
   param->nindicator    =0;
   param->indicator     =NULL;

   // type of multigrid
   if (param->ipar[10]==0)
      param->amg           ="ilu";
   else if (param->ipar[10]==1)
      param->amg           ="amli";
   else if (param->ipar[10]==2)
      param->amg           ="mg";

   param->npresmoothing =param->ipar[12];
   param->npostsmoothing=param->ipar[13];
   param->ncoarse       =param->ipar[11];

   // pre- and post smoother (if used with "mg")
   if (param->ipar[14]==1)
      param->presmoother   ="gsf";
   else if (param->ipar[14]==2)
      param->presmoother   ="gsb";
   else if (param->ipar[14]==3)
      param->presmoother   ="j";
   else if (param->ipar[14]==4)
      param->presmoother   ="ilu";
   if (param->ipar[15]==1)
      param->postsmoother  ="gsf";
   else if (param->ipar[15]==2)
      param->postsmoother  ="gsb";
   else if (param->ipar[15]==3)
      param->postsmoother  ="j";
   else if (param->ipar[15]==4)
      param->postsmoother  ="ilu";

   param->FCpart        ="none";

   // type of coarse grid system
   if (param->ipar[6]&SIMPLE_SC)
      param->typecoarse    ="ilu";
   else if (!(param->ipar[6]&TISMENETSKY_SC))
      param->typecoarse    ="amg";

   param->nrestart      =param->ipar[24];
   param->flags         =param->ipar[6];

   // set up solver and prepare memory allocation for the 
   // iterative solver 
   i=param->ipar[24];
   switch (param->ipar[5]) {
   case  2:  // sbcg
     param->solver="sbcg";
     i=5*A->nr;
     break;
   case  3:  // bcg
     param->solver="bcg";
     i=7*A->nr;
     break;
   case  8:  // gmres 
     param->solver="gmres";
     i=(A->nr+3)*(i+2)+(i+1)*i/2;
     break;
   case  9:  // fGMRES
     param->solver="fgmres";
     i=2*A->nr*(i+1)+((i+1)*i)/2 + 3*i + 2;
     break;
   case 22:  // fsqmr
     param->solver="fsqmr";
     i=(7+2*i)*A->nr;
     break;
   default:  // sqmr 
     param->ipar[5]=4;
     param->solver="sqmr";
     i=6*A->nr;
   } // end switch
   // memory required by the iterative method
   param->ipar[23]=i;

#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
   param->damping=2.0/3.0;
#else
   param->damping=2.0/3.0;
#endif
   param->contraction=0.5;
   param->mixedprecision=0;

   param->nshifts=0;
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
   param->shift0=0;
   param->shiftmax=0;
#else
   param->shift0.r=0;
   param->shift0.i=0;
   param->shiftmax.r=0;
   param->shiftmax.i=0;
#endif
   param->shifts=NULL;
   param->shiftmatrix=NULL;

   param->nthreads=1;
   param->loadbalancefactor=2.0;

   param->indpartial=NULL;
   param->indexpartial=NULL;
   param->valpartial=NULL;
   param->valuepartial=NULL;
   param->vpartial=NULL;
   param->tvpartial=NULL;

   param->nblocks=0;
   param->blocksize=NULL;
} // end symAMGinits



// dummy routine
#else
void myilupackdummy6()
{
}
#endif
