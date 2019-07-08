#!/bin/bash

##########
#name:          set_guix_profile.sh
#description:   load guix profile
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


echo "Loading guix profile"
export GUIX_PROFILE=${HOME}/guix-profiles/framework/.guix-profile
source ${GUIX_PROFILE}/etc/profile

GUIX_LOCPATH="${GUIX_PROFILE}/lib/locale"
unset PYTHONPATH


