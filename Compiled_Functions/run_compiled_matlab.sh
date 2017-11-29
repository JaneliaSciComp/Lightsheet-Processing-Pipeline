#!/bin/bash

ABSOLUTE_SCRIPT=`readlink -m $0`
BIN_DIR=`dirname ${ABSOLUTE_SCRIPT}`

. ${BIN_DIR}/setup_matlab_env.sh

echo """
==============================================================================
  On Host: `hostname`
  Running $*

  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}

  PATH=${PATH}

  XAPPLRESDIR=${XAPPLRESDIR}
  MCR_INHIBIT_CTF_LOCK=${MCR_INHIBIT_CTF_LOCK}
  MCR_CACHE_ROOT=${MCR_CACHE_ROOT}

==============================================================================
"""

$*

if [[ ${MCR_CACHE_ROOT} == /scratch/* ]]; then
  echo """
removing MCR_CACHE_ROOT "${MCR_CACHE_ROOT}" ...
"""
  rm -rf "${MCR_CACHE_ROOT}"
fi
