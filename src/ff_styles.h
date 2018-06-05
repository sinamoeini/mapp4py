#include "ff_lj.h"
#include "ff_eam.h"

#include "ff_fs.h"
#include "ff_eam_dmd.h"
#ifdef SC_DMD
#include "ff_eam_dmd_sc.h"
#include "ff_eam_dmd_scc.h"
#include "ff_eam_dmd_cluster.h"
#endif
#ifdef POTFIT
#include "ff_eam_fit.h"
#include "ff_eam_fit_o.h"
#endif
